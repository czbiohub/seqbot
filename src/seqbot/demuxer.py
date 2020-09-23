#!/usr/bin/env python


# script to:
#   - scan sequencer folders
#   - check for completed runs
#   - demultiplex the run
#   - upload fastq files to S3 when completed

import logging
import os
import pathlib
import subprocess
import sys
import time

import click

import utilities.log_util as ut_log
import utilities.s3_util as s3u

import seqbot.mailbot as mailbot
import seqbot.util as util


log = logging.getLogger(__name__)


config = util.get_config()

# locations where sequencing data gets written
SEQ_DIRS = [pathlib.Path(seq_dir) for seq_dir in config["seqs"]["base"]]

# location on S3 where fastqs are sent
S3_FASTQ_URI = f's3://{config["s3"]["seqbot_bucket"]}/{config["s3"]["fastq_prefix"]}'
S3_SAMPLESHEET_URI = (
    f's3://{config["s3"]["seqbot_bucket"]}/{config["s3"]["samplesheet_prefix"]}'
)

# cache files
index_cache = pathlib.Path(config["local"]["index_cache"])
samplesheet_dir = pathlib.Path(config["local"]["samplesheet_dir"])
output_dir = pathlib.Path(config["local"]["output_dir"])


def maybe_exit_process():
    # get all python pids
    pids = subprocess.run(
        "pgrep demuxer", shell=True, universal_newlines=True, stdout=subprocess.PIPE
    ).stdout.split()

    cmds = 0
    for pid in pids:
        with open(pathlib.Path("/proc") / pid / "cmdline") as f:
            # get the cmdline and match against this script
            line = f.read().split("\x00")[1]
            if line == sys.argv[0]:
                cmds += 1

    # if there are more than one, exit this one
    if cmds > 1:
        print("existing process found, exiting", file=sys.stderr)
        sys.exit(0)
    else:
        # otherwise, take a short nap before starting
        time.sleep(5)


def run_bcl2fastx(seq_dir: pathlib.Path, demux_cmd: list, name: str):
    log.debug(f"running command:\n\t{' '.join(demux_cmd)}")

    try:
        proc = subprocess.run(
            demux_cmd,
            universal_newlines=True,
            capture_output=True,
            timeout=config["demux"]["timeout"],
        )
    except subprocess.TimeoutExpired as exc:
        log.error(f"{name} timed out after {exc.timeout} seconds")
        mailbot.timeout_mail(seq_dir.name, exc.timeout, config["email"], "demux")
        return False

    if proc.returncode != 0:
        log.error(f"{name} returned code {proc.returncode}")
        mailbot.error_mail(seq_dir.name, proc, config["email"], "demux")
        return False
    else:
        return True


def demux_run(seq_dir: pathlib.Path):
    try:
        hdr, rows, split_lanes, cellranger, index_overlap = util.get_samplesheet(
            seq_dir, config, samplesheet_dir / f"{seq_dir.name}.csv"
        )
    except ValueError:
        log.exception("Bad samplesheet:")
        return False
    except UnicodeDecodeError:
        mailbot.samplesheet_error_mail(seq_dir.name, config["email"])
        return False

    if len(rows) == 0:
        log.error("Samplesheet is empty, skipping")
        return False

    samplesheet_path = samplesheet_dir / f"{seq_dir.name}.csv"
    fastq_output = output_dir / seq_dir.name
    fastq_output.mkdir(exist_ok=True)

    log.info(f"demuxing {seq_dir}")

    loading_threads = min(config["demux"]["local_threads"] // 2, len(rows))
    writing_threads = min(config["demux"]["local_threads"] // 2, len(rows))

    demux_cmd = [
        config["demux"]["bcl2fastq"],
        "--processing-threads",
        f"{config['demux']['local_threads']}",
        "--loading-threads",
        f"{loading_threads}",
        "--writing-threads",
        f"{writing_threads}",
        "--sample-sheet",
        f"{samplesheet_path}",
        "--runfolder-dir",
        f"{seq_dir}",
        "--output-dir",
        f"{fastq_output}",
    ]

    if cellranger:
        if cellranger not in (2, 3, "VDJ"):
            log.error("Unknown cellranger version, skipping")
            return False

        demux_cmd.extend(
            [
                "--use-bases-mask=Y*,I*,Y*",
                "--create-fastq-for-index-reads",
                "--minimum-trimmed-read-length=8",
                "--mask-short-adapter-reads=8",
                "--ignore-missing-positions",
                "--ignore-missing-controls",
                "--ignore-missing-filter",
                "--ignore-missing-bcls",
            ]
        )
    elif not split_lanes:
        demux_cmd.append("--no-lane-splitting")

    if index_overlap:
        demux_cmd.extend(["--barcode-mismatches", "0"])

    if len(rows) <= config["demux"]["split_above"]:
        succeeded = run_bcl2fastx(seq_dir, demux_cmd, "bcl2fastq")
        if not succeeded:
            return False
    elif seq_dir.parent.name.startswith("NovaSeq") and not cellranger:
        demux_cmd = [
            config["demux"]["bcl2fastr"],
            "--threads",
            f"{config['demux']['local_threads']}",
            "--read-chunks",
            f"{config['demux']['read_chunks']}",
            "--samplesheet",
            f"{samplesheet_path}",
            "--run-path",
            f"{seq_dir}",
            "--output",
            f"{fastq_output}",
        ]

        succeeded = run_bcl2fastx(seq_dir, demux_cmd, "bcl2fastr")
        if not succeeded:
            return False
    elif not cellranger:
        log.warning(f"{len(rows)} is just too much for now")
        return False
    else:
        log.warning(f"Cellranger run with {len(rows)} samples is not supported")
        return False

    upload_cmd = [
        config["s3"]["awscli"],
        "s3",
        "sync",
        "--no-progress",
        f"{fastq_output}",
        f"{S3_FASTQ_URI}/{seq_dir.name}",
    ]

    log.debug(f"uploading results: '{' '.join(upload_cmd)}'")

    proc = subprocess.run(upload_cmd, universal_newlines=True, capture_output=True)
    failed = proc.returncode != 0

    if failed:
        log.error(f"upload to s3 returned code {proc.returncode}")
        mailbot.error_mail(
            seq_dir.name, proc, email_config=config["email"], task="upload"
        )

        return False
    else:
        if config["local"]["clean"]:
            rm_cmd = ["rm", "-rf", f"{fastq_output}"]
            log.debug(f"removing local copy: '{' '.join(rm_cmd)}'")
            proc = subprocess.run(rm_cmd, universal_newlines=True, capture_output=True)

            if proc.returncode != 0:
                log.warning(f"rm failed for {seq_dir.name}!")
                log.warning(proc.stderr)

        log.info("Sending notification email")
        mailbot.demux_mail(S3_FASTQ_URI, seq_dir.name, config["email"], index_overlap)

        return True


def novaseq_index(seq_dir: pathlib.Path):
    output_file = output_dir / seq_dir.name / f"{seq_dir.name}.txt"
    output_file.parent.mkdir(exist_ok=True)

    index_cmd = [
        config["index"]["bcl2index"],
        "--threads",
        f"{config['index']['local_threads']}",
        "--run-path",
        f"{seq_dir}",
        "--output",
        output_file.as_posix(),
        "--top-n",
        f"{config['index']['index_top_n']}",
    ]

    log.debug(f"running command:\n\t{' '.join(index_cmd)}")

    try:
        proc = subprocess.run(
            index_cmd,
            universal_newlines=True,
            capture_output=True,
            timeout=config["index"]["timeout"],
        )
    except subprocess.TimeoutExpired as exc:
        log.error(f"index timed out after {exc.timeout} seconds")
        mailbot.timeout_mail(seq_dir.name, exc.timeout, config["email"], "index")
        return True

    if proc.returncode != 0:
        log.error(f"bcl2index returned code {proc.returncode}")
        mailbot.error_mail(seq_dir.name, proc, config["email"], "index")
        return False

    log.info("Sending index counts email")
    mailbot.mail_nova_index(seq_dir.name, output_file, config["email"], False)

    if config["local"]["clean"]:
        log.debug(f"Cleaning up {output_file}")
        os.remove(output_file)

    return True


def filter_index(seq_dir: pathlib.Path):
    input_file = output_dir / seq_dir.name / f"{seq_dir.name}.txt"
    output_file = output_dir / seq_dir.name / f"{seq_dir.name}_filtered.txt"
    samplesheet_path = samplesheet_dir / f"{seq_dir.name}.csv"

    if not input_file.exists():
        log.warning(f"{input_file} doesn't exist, skipping")
        return False

    if not samplesheet_path.exists():
        log.warning(f"Can't find {samplesheet_path}, skipping")

    filter_cmd = [
        config["index"]["filter_index"],
        "--input-count",
        input_file.as_posix(),
        "--output",
        output_file.as_posix(),
        "--samplesheet",
        samplesheet_path.as_posix(),
    ]

    log.debug(f"running command:\n\t{' '.join(filter_cmd)}")

    proc = subprocess.run(filter_cmd, universal_newlines=True, capture_output=True)

    if proc.returncode != 0:
        log.error(f"filter_index returned code {proc.returncode}")
        mailbot.error_mail(seq_dir.name, proc, config["email"], "filter")
        return False

    log.info("Sending filtered index counts email")
    mailbot.mail_nova_index(seq_dir.name, output_file, config["email"], True)

    if config["local"]["clean"]:
        log.debug(f"Cleaning up {output_file}")
        os.remove(output_file)

    return True


@click.command()
def main():
    # check for an existing process running
    maybe_exit_process()

    ut_log.get_trfh_logger(
        __package__,
        (config["logging"]["info"], logging.INFO, "midnight", 14),
        (config["logging"]["debug"], logging.DEBUG, "midnight", 3),
    )

    log.debug("querying s3 for list of existing runs")
    demux_set = {
        os.path.basename(d[:-1])
        for d in s3u.get_folders(
            bucket=config["s3"]["seqbot_bucket"],
            prefix=config["s3"]["fastq_prefix"] + "/",
        )
    }

    if index_cache.exists():
        log.debug("reading cache file for indexed runs")
        with index_cache.open() as f:
            index_set = {line.strip() for line in f}
    else:
        index_set = set()

    log.info(f"{len(demux_set)} folders in {S3_FASTQ_URI}")

    log.debug("Getting the list of sample-sheets")
    samplesheets = {
        os.path.splitext(os.path.basename(fn))[0]
        for fn in s3u.get_files(
            bucket=config["s3"]["seqbot_bucket"],
            prefix=config["s3"]["samplesheet_prefix"],
        )
    }

    log.info(f"{len(samplesheets)} samplesheets in {S3_SAMPLESHEET_URI}")

    updated_demux_set = demux_set.copy()
    updated_index_set = index_set.copy()

    # for each sequencer, check for newly completed runs
    for seq_dir in SEQ_DIRS:
        log.info(f"scanning {seq_dir}...")

        for seq in config["seqs"]["dirs"]:
            log.info(seq)
            fns = sorted(
                (seq_dir / seq).glob(f'[0-9]*/{config["seqs"]["sentinels"][seq]}')
            )
            log.debug(f"{len(fns)} ~complete runs in {seq_dir / seq}")

            for fn in fns:
                run_dir = fn.parent
                if run_dir.name in updated_demux_set:
                    continue

                if not seq.startswith("NovaSeq"):
                    if run_dir.name not in samplesheets:
                        log.debug(f"skipping {run_dir.name}, no sample-sheet")
                        continue

                    if demux_run(run_dir):
                        updated_demux_set.add(run_dir.name)

                # elif (run_dir / "CopyComplete.txt").exists():
                else:
                    if run_dir.name not in index_set:
                        log.debug(f"no record for {run_dir.name}, indexing...")
                        if novaseq_index(run_dir):
                            updated_index_set.add(run_dir.name)

                    if run_dir.name not in samplesheets:
                        log.debug(f"skipping {run_dir.name}, no sample-sheet")
                        continue
                    elif demux_run(run_dir):
                        updated_demux_set.add(run_dir.name)

                        if run_dir.name in updated_index_set:
                            log.debug(f"filtering index count for {run_dir.name}")
                            filter_index(run_dir)

    log.info("scan complete")

    log.info(f"demuxed {len(updated_demux_set) - len(demux_set)} new runs")
    log.info(f"indexed {len(updated_index_set) - len(index_set)} new runs")

    with index_cache.open(mode="w") as OUT:
        print("\n".join(sorted(updated_index_set)), file=OUT)

    log.debug("wrote new cache file")
