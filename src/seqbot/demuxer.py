#!/usr/bin/env python


# script to:
#   - scan sequencer folders
#   - check for completed runs
#   - demultiplex the run
#   - upload fastq files to S3 when completed


import collections
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
scratch_space = pathlib.Path(config["local"]["scratch"])


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


def run_bcl2fastq(seq_dir: pathlib.Path, demux_cmd: list, logger: logging.Logger):
    logger.debug(f"running command:\n\t{' '.join(demux_cmd)}")

    try:
        proc = subprocess.run(
            demux_cmd,
            universal_newlines=True,
            capture_output=True,
            timeout=config["demux"]["timeout"]
        )
    except subprocess.TimeoutExpired as exc:
        logger.error(f"bcl2fastq timed out after {exc.timeout} seconds")
        logger.debug(f"sending mail to:")
        logger.debug(",".join(config["email"]["addresses_to_email_on_error"]))
        mailbot.timeout_mail(seq_dir.name, exc.timeout, config["email"])
        return False

    if proc.returncode != 0:
        logger.error(f"bcl2fastq returned code {proc.returncode}")
        logger.debug(f"sending error mail to:")
        logger.debug(",".join(config["email"]["addresses_to_email_on_error"]))
        mailbot.error_mail(seq_dir.name, proc, config["email"])
        return False
    else:
        return True


def demux_run(seq_dir: pathlib.Path, logger: logging.Logger):
    try:
        hdr, rows, split_lanes, cellranger, index_overlap = util.get_samplesheet(
            seq_dir, config, logger, samplesheet_dir / f"{seq_dir.name}.csv"
        )
    except ValueError:
        return False
    except UnicodeDecodeError:
        logger.debug("Sending error mail to:")
        logger.debug(",".join(config["email"]["addresses_to_email"]))
        mailbot.samplesheet_error_mail(seq_dir.name, config["email"])
        return False

    samplesheet_path = samplesheet_dir / f"{seq_dir.name}.csv"
    temp_output = scratch_space / seq_dir.name
    temp_output.mkdir(exist_ok=True)

    logger.info(f"demuxing {seq_dir}")

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
        f"{temp_output}",
    ]

    if cellranger:
        if cellranger == 2:
            demux_cmd.append("--use-bases-mask=Y26,I8,Y98")
        elif cellranger == 3:
            demux_cmd.append("--use-bases-mask=Y28,I8,Y91")
        else:
            logger.error("Unknown cellranger version, skipping")
            return False

        demux_cmd.extend(
            [
                "--create-fastq-for-index-reads",
                "--minimum-trimmed-read-length=8",
                "--mask-short-adapter-reads=8",
                "--ignore-missing-positions",
                "--ignore-missing-controls",
                "--ignore-missing-filter",
                "--ignore-missing-bcls",
            ]
        )


    if not split_lanes:
        demux_cmd.append("--no-lane-splitting")

    if index_overlap:
        demux_cmd.extend(["--barcode-mismatches", "0"])

    if len(rows) <= config["demux"]["split_above"]:
        succeeded = run_bcl2fastq(seq_dir, demux_cmd, logger)
        if not succeeded:
            return False
    else:
        logger.info("Not demuxing big runs yet")
        return False

    upload_cmd = [
        config["s3"]["awscli"],
        "s3",
        "sync",
        "--no-progress",
        f"{temp_output}",
        f"{S3_FASTQ_URI}/{seq_dir.name}",
    ]

    logger.debug(f"uploading results: '{' '.join(upload_cmd)}'")

    proc = subprocess.run(upload_cmd, universal_newlines=True, capture_output=True)
    failed = proc.returncode != 0

    if failed:
        logger.error(f"upload to s3 returned code {proc.returncode}")
        logger.debug(f"sending error mail to:")
        logger.debug(",".join(config["email"]["addresses_to_email_on_error"]))
        mailbot.error_mail(seq_dir.name, proc, email_config=config["email"])

        return False
    else:
        if config["local"]["clean"]:
            rm_cmd = ["rm", "-rf", f"{temp_output}"]
            logger.debug(f"removing local copy: '{' '.join(rm_cmd)}'")
            proc = subprocess.run(rm_cmd, universal_newlines=True, capture_output=True)

            if proc.returncode != 0:
                logger.warning(f"rm failed for {seq_dir.name}!")
                logger.warning(proc.stderr)

        logger.info("Sending notification email")
        mailbot.demux_mail(S3_FASTQ_URI, seq_dir.name, config["email"], index_overlap)

        return True


def novaseq_index(seq_dir: pathlib.Path, logger: logging.Logger):
    output_file = scratch_space / f"{seq_dir.name}.txt"

    index_cmd = [
        config["demux"]["nova_index"],
        "-p",
        f"{config['demux']['local_threads']}",
        "--run_path",
        f"{seq_dir}",
        "--output_file",
        output_file.as_posix(),
        "--top_n_counts",
        f"{config['demux']['index_top_n']}",
    ]

    logger.debug(f"running command:\n\t{' '.join(index_cmd)}")

    proc = subprocess.run(index_cmd, universal_newlines=True, capture_output=True)

    if proc.returncode != 0:
        logger.error(f"nova_index returned code {proc.returncode}")
        logger.debug(f"sending error mail to:")
        logger.debug(",".join(config["email"]["addresses_to_email_on_error"]))
        mailbot.error_mail(seq_dir.name, proc, config["email"])
        return False

    logger.info("Sending index counts email")
    mailbot.mail_nova_index(seq_dir.name, output_file, config["email"])

    if config["local"]["clean"]:
        logger.debug(f"Cleaning up {output_file}")
        os.remove(output_file)

    return True


@click.command()
def main():
    # check for an existing process running
    maybe_exit_process()

    logger = ut_log.get_trfh_logger(
        "demuxer",
        (config["logging"]["info"], logging.INFO, "midnight", 14),
        (config["logging"]["debug"], logging.DEBUG, "midnight", 3),
    )

    logger.debug("querying s3 for list of existing runs")
    demux_set = {
        os.path.basename(d[:-1])
        for d in s3u.get_folders(
            bucket=config["s3"]["seqbot_bucket"],
            prefix=config["s3"]["fastq_prefix"] + "/",
        )
    }

    if index_cache.exists():
        logger.debug("reading cache file for indexed runs")
        with index_cache.open() as f:
            index_set = {line.strip() for line in f}
    else:
        index_set = set()

    logger.info(f"{len(demux_set)} folders in {S3_FASTQ_URI}")

    logger.debug("Getting the list of sample-sheets")
    samplesheets = {
        os.path.splitext(os.path.basename(fn))[0]
        for fn in s3u.get_files(
            bucket=config["s3"]["seqbot_bucket"],
            prefix=config["s3"]["samplesheet_prefix"],
        )
    }

    logger.info(f"{len(samplesheets)} samplesheets in {S3_SAMPLESHEET_URI}")

    updated_demux_set = demux_set.copy()
    updated_index_set = index_set.copy()

    # for each sequencer, check for newly completed runs
    for seq_dir in SEQ_DIRS:
        logger.info(f"scanning {seq_dir}...")

        for seq in config["seqs"]["dirs"]:
            logger.info(seq)
            fns = sorted(
                (seq_dir / seq).glob(f'[0-9]*/{config["seqs"]["sentinels"][seq]}')
            )
            logger.debug(f"{len(fns)} ~complete runs in {seq_dir / seq}")

            for fn in fns:
                run_dir = fn.parent
                if run_dir.name in updated_demux_set:
                    logger.debug(f"skipping {run_dir.name}, already demuxed")
                    continue

                if not seq.startswith("NovaSeq"):
                    if run_dir.name not in samplesheets:
                        logger.debug(f"skipping {run_dir.name}, no sample-sheet")
                        continue

                    if demux_run(run_dir, logger):
                        updated_demux_set.add(run_dir.name)

                elif (run_dir / "CopyComplete.txt").exists():
                    if (
                        run_dir.name not in samplesheets
                        and run_dir.name not in index_set
                    ):
                        logger.debug(f"no samplesheet for {run_dir.name}, indexing...")
                        if novaseq_index(run_dir, logger):
                            updated_index_set.add(run_dir.name)
                    elif run_dir.name not in samplesheets:
                        logger.debug(f"skipping {run_dir.name}, no sample-sheet")
                        continue
                    elif demux_run(run_dir, logger):
                        updated_demux_set.add(run_dir.name)

    logger.info("scan complete")

    logger.info(f"demuxed {len(updated_demux_set) - len(demux_set)} new runs")
    logger.info(f"indexed {len(updated_index_set) - len(index_set)} new runs")

    with index_cache.open(mode="w") as OUT:
        print("\n".join(sorted(updated_index_set)), file=OUT)

    logger.debug("wrote new cache file")
