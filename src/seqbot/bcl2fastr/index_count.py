#!/usr/bin/env python

import gzip
import itertools
import logging
import math
import pathlib
import traceback

from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing import Manager, Queue

import click

import numpy as np

from seqbot import bcl2fastr

import utilities.log_util as ut_log


def index_count_thread(lane: int, part: int, i: int):
    read_to_indexes = lambda r: "\t".join(
        r[i0:i1] for i0, i1 in zip(index_splits, index_splits[1:])
    )

    read_counter = Counter(
        read_to_indexes(read)
        for read, qscore, read_id in bcl2fastr.extract_reads(novaseq_run, lane, part, i)
    )

    return read_counter


def index_count_processor(
    i: int, n_proc: int, n_threads: int, top_n: int, log_queue: Queue
):
    try:
        log_queue.put(
            (
                f"starting job with args: ({i}, {n_proc}, {n_threads}, {top_n})",
                logging.DEBUG,
            )
        )

        tile_gen = (
            (lane, part, i)
            for lane, part in novaseq_run.lane_parts
            for i in range(novaseq_run.tiles[lane, part].shape[0])
        )

        total_counts = Counter()

        log_queue.put((f"Starting thread pool of {n_threads} threads", logging.DEBUG))
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            for top_counts in executor.map(
                index_count_thread, *zip(*itertools.islice(tile_gen, i, None, n_proc))
            ):
                total_counts.update(top_counts)

        log_queue.put(
            (f"job done for args: ({i}, {n_proc}, {n_threads}, {top_n})", logging.DEBUG)
        )

        return Counter(dict(total_counts.most_common(top_n)))
    except Exception as exc:
        log_queue.put(("".join(traceback.format_tb(exc.__traceback__)), logging.ERROR))
        log_queue.put(
            ("encountered exception in process:\n{}".format(exc), logging.ERROR)
        )
        return tuple()


@click.command()
@click.option("-p", "--n_processes", required=True, type=int)
@click.option(
    "-r",
    "--run_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option("-o", "--output_file", required=True, type=click.Path())
@click.option("-t", "--n_threads", type=int, default=8)
@click.option("--top_n_counts", type=int, default=384)
@click.option("-d", "--debug", is_flag=True)
def main(
    n_processes: int,
    run_path: str,
    output_file: str,
    n_threads: int,
    top_n_counts: int,
    debug: bool,
):
    logger, _, _ = ut_log.get_logger("index_count", debug=debug)

    run_path = pathlib.Path(run_path)
    output_file = pathlib.Path(output_file)

    logger.info(f"Reading headers and filters for {run_path}")

    run_info = bcl2fastr.NovaSeqRun.read_run_info(run_path)

    index_reads = [(i0, i1, is_index) for i0, i1, is_index in run_info if is_index]

    global novaseq_run
    novaseq_run = bcl2fastr.NovaSeqRun(
        run_path, start=index_reads[0][0] + 1, stop=index_reads[-1][1] + 1
    )

    global index_splits
    index_splits = (0,) + tuple(
        np.cumsum([(i1 - i0) for i0, i1, is_index in index_reads])
    )

    # each process job will take enough jobs to populate its threads once
    n_tiles = sum(
        novaseq_run.tiles[lane, part].shape[0] for lane, part in novaseq_run.lane_parts
    )
    n_proc_jobs = math.ceil(n_tiles / n_threads)

    manager = Manager()

    log_queue, log_thread = ut_log.get_thread_logger(logger, queue=manager.Queue())

    n_files = sum(map(len, novaseq_run.cbcl_files.values()))

    logger.info(f"reading {n_files} files and aggregating counters")

    total_counts = Counter()

    logger.debug(f"initializing pool of {n_processes} processes")
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        logger.debug("starting count")

        for top_counts in executor.map(
            index_count_processor,
            range(n_proc_jobs),
            itertools.repeat(n_proc_jobs),
            itertools.repeat(n_threads),
            itertools.repeat(8 * top_n_counts),
            itertools.repeat(log_queue),
        ):
            total_counts.update(top_counts)

    log_queue.put("STOP")
    log_thread.join()

    if output_file.suffix == ".gz":
        with gzip.open(output_file, "wt") as out:
            for index, count in total_counts.most_common(top_n_counts):
                print(f"{index}\t{count}", file=out)
    else:
        with open(output_file, "w") as out:
            for index, count in total_counts.most_common(top_n_counts):
                print(f"{index}\t{count}", file=out)

    logger.info("done!")
