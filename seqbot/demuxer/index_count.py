#!/usr/bin/env python

import gzip
import itertools
import logging
import pathlib

from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor

import click

import numpy as np

import seqbot.demuxer.bcl2fastr as bcl2fastr

import utilities.log_util as ut_log


def index_count_processor(
    cbcl_files, cbcl_filter_files, loc_file, i, nproc, top_n, i_splits
):
    read_to_indexes = lambda r: "\t".join(
        r[i0:i1] for i0, i1 in zip(i_splits, i_splits[1:])
    )

    try:
        msg = (
            f"starting job with args: "
            f"({cbcl_files[0]}..., {cbcl_filter_files[0]}..., {loc_file}, {i}, {nproc})"
        )
        log_queue.put((msg, logging.DEBUG))

        read_counter = Counter(
            read_to_indexes(read)
            for read, qscore, read_id in bcl2fastr.extract_reads(
                cbcl_files, cbcl_filter_files, loc_file, i, nproc
            )
        )

        msg = (
            f"job done for args: "
            f"({cbcl_files[0]}..., {cbcl_filter_files[0]}..., {loc_file}, {i}, {nproc})"
        )
        log_queue.put((msg, logging.DEBUG))

        return Counter(dict(read_counter.most_common(top_n)))
    except Exception as detail:
        log_queue.put(
            ("encountered exception in process:\n{}".format(detail), logging.INFO)
        )
        return tuple()


@click.command()
@click.option("--n_threads", required=True, type=int)
@click.option(
    "--bcl_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option("--output_file", required=True, type=click.Path())
@click.option("--debug", type=bool)
@click.option("--top_n_counts", type=int)
def main(n_threads, bcl_path, output_file, debug=False, top_n_counts=20000):
    logger, _, _ = ut_log.get_logger("index_count", debug=debug)

    bcl_path = pathlib.Path(bcl_path)

    cbcl_file_lists, cbcl_filter_lists, loc_file = bcl2fastr.cbcl_globber(bcl_path)

    read_info = bcl2fastr.get_read_info(bcl_path / "RunInfo.xml")

    index_set = {j for i0, i1, is_index in read_info if is_index for j in range(i0, i1)}
    index_splits = (0,) + tuple(
        np.cumsum([(i1 - i0) for i0, i1, is_index in read_info if is_index])
    )

    in_range = lambda cfn: bcl2fastr.get_cycle(cfn) in index_set

    cbcl_file_lists = {
        (lane, part): tuple(cfn for cfn in cbcl_file_lists[lane, part] if in_range(cfn))
        for lane, part in cbcl_file_lists
    }

    logger.info("{} CBCL files to read".format(sum(map(len, cbcl_file_lists.values()))))

    lane_parts = sorted(cbcl_file_lists)

    global log_queue
    log_queue, log_thread = ut_log.get_thread_logger(logger)

    logger.info(
        "reading {} files and aggregating counters".format(
            sum(map(len, cbcl_file_lists.values()))
        )
    )

    # warning: gratuitous use of itertools module ahead! it's gonna be great

    # lambda function to make this crazy itertools chain.
    # looks nuts, just repeats each element of s for [n_threads] times
    rep_n = lambda s: itertools.chain.from_iterable(
        map(itertools.repeat, s, itertools.repeat(n_threads))
    )

    total_counts = Counter()

    logger.debug(f"initializing pool of {n_threads} processes")
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        logger.debug("starting count")
        for i, top_counts in enumerate(
            executor.map(
                index_count_processor,
                rep_n(cbcl_file_lists[lane, part] for lane, part in lane_parts),
                rep_n(cbcl_filter_lists[lane] for lane, part in lane_parts),
                itertools.repeat(loc_file),
                itertools.cycle(range(n_threads)),
                itertools.repeat(n_threads),
                itertools.repeat(8 * top_n_counts),
                itertools.repeat(index_splits),
            )
        ):
            if i % 100 == 0:
                logger.info(f"{i}")
            total_counts.update(top_counts)

    log_queue.put("STOP")
    log_thread.join()

    with gzip.open(output_file, "wt") as out:
        for index, count in total_counts.most_common(top_n_counts):
            print(f"{index}\t{count}", file=out)

    logger.info("done!")
