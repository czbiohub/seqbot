#!/usr/bin/env python

import gzip
import itertools
import logging
import pathlib

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import click

import numpy as np

import seqbot.demuxer.bcl2fastr as bcl2fastr

import utilities.log_util as ut_log


def read_processor(
    cbcl_files, cbcl_filter_files, loc_file, i, nproc, read_tmp, out_file
):
    try:
        msg = (
            f"starting job with args: "
            f"({cbcl_files[0]}..., {cbcl_filter_files[0]}..., {loc_file}, {i}, {nproc})"
        )
        log_queue.put((msg, logging.DEBUG))

        with gzip.open(out_file, "wt") as OUT:
            for read, qscore, read_id in bcl2fastr.extract_reads(
                cbcl_files, cbcl_filter_files, loc_file, i, nproc
            ):
                print(read_tmp.format(read_id), file=OUT)
                print(read, file=OUT)
                print("+", file=OUT)
                print(qscore, file=OUT)

        msg = (
            f"job done for args: "
            f"({cbcl_files[0]}..., {cbcl_filter_files[0]}..., {loc_file}, {i}, {nproc})"
        )
        log_queue.put((msg, logging.DEBUG))
    except Exception as detail:
        log_queue.put(
            ("encountered exception in process:\n{}".format(detail), logging.INFO)
        )


@click.command()
@click.option("--n_threads", required=True, type=int)
@click.option(
    "--bcl_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option(
    "--samplesheet",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--output_dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option("--debug", type=bool)
def main(n_threads, bcl_path, samplesheet, output_dir, debug=False):
    logger, _, _ = ut_log.get_logger("write_fastq", debug=debug)

    bcl_path = pathlib.Path(bcl_path)
    samplesheet = pathlib.Path(samplesheet)

    cbcl_file_lists, cbcl_filter_lists, loc_file = bcl2fastr.cbcl_globber(bcl_path)

    read_info = bcl2fastr.get_read_info(bcl_path / "RunInfo.xml")

    in_range = lambda cfn: (
        index_cycle_start <= bcl2fastr.get_cycle(cfn) < index_cycle_end
    )

    cbcl_file_lists = {
        (lane, part): tuple(cfn for cfn in cbcl_file_lists[lane, part] if in_range(cfn))
        for lane, part in cbcl_file_lists
    }

    logger.info("{} CBCL files to read".format(sum(map(len, cbcl_file_lists.values()))))

    lane_parts = sorted(cbcl_file_lists)

    global log_queue
    log_queue, log_thread = ut_log.get_thread_logger(logger)

    logger.info(
        "reading {} files and writing fastq.gz files".format(
            sum(map(len, cbcl_file_lists.values()))
        )
    )

    rid_tmp = f"{bcl2fastr.cbcl_id(bcl_path)}:{{}}:{{}} 1:N:0:0"
    output_file = str(output_dir / "read_file_{}.fastq.gz")

    # warning: gratuitous use of itertools module ahead! it's gonna be great

    # lambda function to make this crazy itertools chain.
    # looks nuts, just repeats each element of s for [n_threads] times
    rep_n = lambda s: itertools.chain.from_iterable(
        map(itertools.repeat, s, itertools.repeat(n_threads))
    )

    logger.debug("initializing pool of {} processes".format(n_threads))
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        logger.debug("starting demux")
        for i, _ in enumerate(
            executor.map(
                read_processor,
                rep_n(cbcl_file_lists[lane, part] for lane, part in lane_parts),
                rep_n(cbcl_filter_lists[lane] for lane, part in lane_parts),
                itertools.repeat(loc_file),
                itertools.cycle(range(n_threads)),
                itertools.repeat(n_threads),
                rep_n(rid_tmp.format(lane, "{}") for lane, part in lane_parts),
                map(output_file.format, itertools.count()),
            )
        ):
            if i % 100 == 0:
                logger.info(f"{i}")

    log_queue.put("STOP")
    log_thread.join()

    logger.info("done!")
