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

from logging.handlers import TimedRotatingFileHandler

import click

import utilities.log_util as ut_log
import utilities.s3_util as s3u

import seqbot.mailbot as mailbot
import seqbot.util as util


config = util.get_config()

# location where sequencer data gets written
SEQ_DIR = pathlib.Path(config["seqs"]["base"])

# location on S3 where fastqs are sent
S3_URI = f's3://{config["s3"]["seqbot_bucket"]}/{config["s3"]["fastq_prefix"]}'

# number of samples to batch when dealing with very large runs
sample_n = config["demux"]["sample_n"]

# cache files
samplesheet_dir = pathlib.Path(config["local"]["samplesheet_dir"])
demux_cache = pathlib.Path(config["local"]["demux_cache"])
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


def demux_run(seq_dir: pathlib.Path, logger: logging.Logger):
    try:
        hdr, rows, split_lanes, cellranger = util.get_samplesheet(
            seq_dir, config, logger
        )
    except ValueError:
        return False

    batched = len(rows) > sample_n

    if cellranger and (split_lanes or batched):
        logger.warning("Cellranger workflow won't use extra options")

    batch_range = range(0, len(rows) + int(len(rows) % sample_n > 0), sample_n)

    for i, j in enumerate(batch_range):
        with open(samplesheet_dir / f"{seq_dir.name}_{i}.csv", "w") as OUT:
            print(hdr, file=OUT)
            for r in rows[j : j + sample_n]:
                print(",".join(r), file=OUT)

    for i in range(len(batch_range)):
        logger.info(f"demuxing batch {i} of {seq_dir}")
        demux_cmd = [
            config["demux"]["reflow"],
            "run",
            "-local",
            "-localdir",
            scratch_space.as_posix(),
            config["demux"]["reflow_demux"],
            "-bcl_path",
            f"local{seq_dir.as_uri()}",
            "-output_path",
            f"{S3_URI}/{seq_dir.name}",
            "-samplesheet",
            f"local{(samplesheet_dir / seq_dir.name).as_uri()}_{i}.csv",
        ]

        if not split_lanes:
            demux_cmd.append("-no_lane_splitting")

        if batched:
            demux_cmd.extend(("-batch_runID", str(i + 1)))

        if cellranger:
            demux_cmd.append("-cellranger")

        logger.debug(f"running command:\n\t{' '.join(demux_cmd)}")

        proc = subprocess.run(
            " ".join(demux_cmd),
            shell=True,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        if proc.returncode != 0:
            logger.error(f"reflow batch {i} returned code {proc.returncode}")
            logger.debug(f"sending error mail to:")
            logger.debug(",".join(config["email"]["addresses_to_email_on_error"]))
            mailbot.error_mail(seq_dir.name, proc, config["email"])
            break
    else:
        logger.info("Sending notification email")
        mailbot.demux_mail(S3_URI, seq_dir.name, config["email"])
        return True

    return False


def novaseq_index(seq_dir: pathlib.Path, logger: logging.Logger):
    output_file = scratch_space / f"{seq_dir.name}.txt.gz"

    index_cmd = [
        config["demux"]["nova_index"],
        "--n_threads",
        f"{config['demux']['local_threads']}",
        "--bcl_path",
        f"{seq_dir}",
        "--output_file",
        output_file.as_posix(),
        "--top_n_counts",
        f"{config['demux']['index_top_n']}",
    ]

    logger.debug(f"running command:\n\t{' '.join(index_cmd)}")

    proc = subprocess.run(
        " ".join(index_cmd),
        shell=True,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    if proc.returncode != 0:
        logger.error(f"nova_index returned code {proc.returncode}")
        logger.debug(f"sending error mail to:")
        logger.debug(",".join(config["email"]["addresses_to_email_on_error"]))
        mailbot.error_mail(seq_dir.name, proc, config["email"])
        return

    logger.info("Sending index counts email")
    mailbot.mail_nova_index(seq_dir.name, output_file, config["email"])

    logger.debug(f"Cleaning up {output_file}")
    os.remove(output_file)


def demux_novaseq(seq_dir: pathlib.Path, logger: logging.Logger):
    return False


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
            bucket=config["s3"]["seqbot_bucket"], prefix=config["s3"]["fastq_prefix"]
        )
    }

    logger.info(f"{len(demux_set)} folders in {S3_URI}")

    logger.debug("Getting the list of sample-sheets")
    samplesheets = {
        os.path.splitext(os.path.basename(fn))[0]
        for fn in s3u.get_files(
            bucket=config["s3"]["seqbot_bucket"], prefix="sample-sheets"
        )
    }

    logger.info(
        f"{len(samplesheets)} samplesheets in "
        f's3://{config["s3"]["seqbot_bucket"]}/sample-sheets'
    )

    updated_demux_set = demux_set.copy()

    logger.info("scanning {}...".format(SEQ_DIR))

    # for each sequencer, check for newly completed runs
    for seq in config["seqs"]["dirs"]:
        logger.info(seq)
        fns = list((SEQ_DIR / seq).glob(f'[0-9]*/{config["seqs"]["sentinels"][seq]}'))
        logger.debug(f"{len(fns)} ~complete runs in {SEQ_DIR / seq}")

        for fn in fns:
            seq_dir = fn.parent
            if seq_dir.name in updated_demux_set:
                logger.debug(f"skipping {seq_dir.name}, already demuxed")
                continue

            if seq != "NovaSeq-01":
                if seq_dir.name not in samplesheets:
                    logger.debug(f"skipping {seq_dir.name}, no sample-sheet")
                    continue

                if demux_run(seq_dir, logger):
                    updated_demux_set.add(seq_dir.name)

            elif (seq_dir / "CopyComplete.txt").exists():
                if seq_dir.name not in samplesheets:
                    logger.debug(f"no samplesheet for {seq_dir.name}, indexing...")
                    novaseq_index(seq_dir, logger)
                    # going to add for now...even though that's not true
                    updated_demux_set.add(seq_dir.name)
                elif demux_novaseq(seq_dir, logger):
                    updated_demux_set.add(seq_dir.name)

    logger.info("scan complete")

    logger.info("demuxed {} new runs".format(len(updated_demux_set) - len(demux_set)))
