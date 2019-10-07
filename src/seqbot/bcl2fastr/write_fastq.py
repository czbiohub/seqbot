#!/usr/bin/env python

import gzip
import itertools
import logging
import math
import os
import pathlib
import traceback

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing import Manager, Queue
from typing import Optional, Sequence, Tuple

import click

import networkx as nx

from seqbot import bcl2fastr
from seqbot.util import hamming_set, read_samplesheet

import utilities.log_util as ut_log


Reads = Tuple[Tuple[str], Tuple[str], Tuple[str]]
Sample = Tuple[str, Optional[int]]
SampleIndex = Tuple[Sample, str]


class SampleData(object):
    def __init__(
        self,
        samplesheet: pathlib.Path,
        output_dir: pathlib.Path,
        novaseq_run: bcl2fastr.NovaSeqRun,
    ):
        with samplesheet.open() as f:
            hdr, h_row, rows = read_samplesheet(f)

        h_row_d = {c: i for i, c in enumerate(h_row)}

        sample_name = h_row_d["sample_name"]
        lane = h_row_d.get("lane", -1)
        self.split_lanes = lane > -1

        self.samples = set()
        self.index_graph = defaultdict(nx.Graph)

        run_lanes = {lane for lane, part in novaseq_run.lane_parts}

        for r in rows:
            sample_lane = int(r[lane]) if self.split_lanes else None
            if sample_lane not in run_lanes:
                continue

            sample = (r[sample_name], sample_lane)
            if "index2" in h_row_d:
                sample_idx = (sample, f"{r[h_row_d['index']]}+{r[h_row_d['index2']]}")
            else:
                sample_idx = (sample, r[h_row_d["index"]])

            self.samples.add(sample)

            self.index_graph[sample_lane].add_node(sample_idx, bipartite=0)
            for i1 in hamming_set(r[h_row_d["index"]]):
                self.index_graph[sample_lane].add_node((i1, 1), bipartite=1)
                self.index_graph[sample_lane].add_edge(sample_idx, (i1, 1))
            if "index2" in h_row_d:
                for i2 in hamming_set(r[h_row_d["index2"]]):
                    self.index_graph[sample_lane].add_node((i2, 2), bipartite=1)
                    self.index_graph[sample_lane].add_edge(sample_idx, (i2, 2))

        if self.split_lanes:
            self.output_files = {
                (sample_name, lane): [
                    (output_dir / f"{sample_name}_L{lane:03}_R{i + 1}_001.fastq.gz")
                    for i in novaseq_run.read_n
                ]
                for sample_name, lane in self.samples
            }
        else:
            self.output_files = {
                (sample_name, lane): [
                    (output_dir / f"{sample_name}_R{i + 1}_001.fastq.gz")
                    for i in novaseq_run.read_n
                ]
                for sample_name, lane in self.samples
            }

    def index_to_sample(
        self, index1: str, index2: Optional[str] = None, lane: Optional[int] = None
    ):
        if (index1, 1) in self.index_graph[lane]:
            if index2 is None:
                sample_set = set(self.index_graph[lane][(index1, 1)])
            elif (index2, 2) in self.index_graph[lane]:
                sample_set = set(self.index_graph[lane][(index1, 1)]) & set(
                    self.index_graph[lane][(index2, 2)]
                )
            else:
                sample_set = set()
        else:
            sample_set = set()

        if len(sample_set) == 1:
            return sample_set.pop()
        else:
            return None


def reader_thread(lane: int, part: int, i: int):
    reads = defaultdict(list)

    for read, qscore, read_id in bcl2fastr.extract_reads(novaseq_run, lane, part, i):
        sample_idx = sample_data.index_to_sample(*novaseq_run.get_indexes(read), lane)

        if sample_idx is not None:
            sample, indexes = sample_idx

            reads[sample].append(
                (
                    novaseq_run.get_reads(read),
                    novaseq_run.get_reads(qscore),
                    tuple(
                        f"{read_id} {i + 1}:N:0:{indexes}" for i in novaseq_run.read_n
                    ),
                )
            )

    return reads


def reader_process(i: int, n_proc: int, n_threads: int, log_queue: Queue):
    try:
        tile_gen = (
            (lane, part, i)
            for lane, part in novaseq_run.lane_parts
            for i in range(novaseq_run.tiles[lane, part].shape[0])
        )

        process_reads = defaultdict(list)

        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            for thread_reads in executor.map(
                reader_thread, *zip(*itertools.islice(tile_gen, i, None, n_proc))
            ):
                for sample, reads in thread_reads.items():
                    process_reads[sample].extend(reads)

        return process_reads
    except Exception as exc:
        log_queue.put(("".join(traceback.format_tb(exc.__traceback__)), logging.ERROR))
        log_queue.put(
            ("encountered exception in process:\n{}".format(exc), logging.ERROR)
        )


def writer_process(sample: Sample, sample_reads: Sequence[Reads], log_queue: Queue):
    try:
        output_files = sample_data.output_files[sample]

        for i, reads in enumerate(zip(*(tuple(zip(*rqr)) for rqr in sample_reads))):
            with gzip.open(output_files[i], "at") as OUT:
                for read, qscore, read_id in reads:
                    print(read_id, read, "+", qscore, sep="\n", file=OUT)
    except Exception as exc:
        log_queue.put(("".join(traceback.format_tb(exc.__traceback__)), logging.ERROR))
        log_queue.put(
            ("encountered exception in process:\n{}".format(exc), logging.ERROR)
        )


@click.command()
@click.option("-p", "--n_processes", required=True, type=int)
@click.option(
    "-r",
    "--run_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option(
    "-s",
    "--samplesheet",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "-o",
    "--output_dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option("-t", "--n_threads", type=int, default=8)
@click.option("-d", "--debug", is_flag=True)
def main(
    n_processes: int,
    run_path: str,
    samplesheet: str,
    output_dir: str,
    n_threads: int,
    debug: bool,
):
    logger, _, _ = ut_log.get_logger("write_fastq", debug=debug)

    run_path = pathlib.Path(run_path)
    samplesheet = pathlib.Path(samplesheet)
    output_dir = pathlib.Path(output_dir)

    logger.info(f"Reading headers and filters for {run_path}")

    global novaseq_run
    novaseq_run = bcl2fastr.NovaSeqRun(run_path)

    global sample_data
    sample_data = SampleData(samplesheet, output_dir, novaseq_run)

    # subset to lanes in the samplesheet
    lanes = {lane for sample_name, lane in sample_data.samples}
    if lanes != {None}:
        novaseq_run.lane_parts = [
            (lane, part) for lane, part in novaseq_run.lane_parts if lane in lanes
        ]

    # each process job will take enough jobs to populate its threads once
    n_tiles = sum(
        novaseq_run.tiles[lane, part].shape[0] for lane, part in novaseq_run.lane_parts
    )
    n_proc_jobs = math.ceil(n_tiles / n_threads)

    manager = Manager()

    log_queue, log_thread = ut_log.get_thread_logger(logger, queue=manager.Queue())

    n_files = sum(map(len, novaseq_run.cbcl_files.values()))

    logger.info(f"reading {n_files} files and writing fastq.gz files")

    for fn in output_dir.glob("*.fastq.gz"):
        logger.warn(f"Removing existing file {fn}")
        os.remove(fn)

    logger.debug("initializing pool of {} processes".format(n_processes))
    with ProcessPoolExecutor(max_workers=n_processes) as read_executor:
        logger.debug("starting demux")

        with ProcessPoolExecutor(max_workers=n_processes) as write_executor:

            for proc_reads in read_executor.map(
                reader_process,
                range(n_proc_jobs),
                itertools.repeat(n_proc_jobs),
                itertools.repeat(n_threads),
                itertools.repeat(log_queue),
            ):
                list(
                    write_executor.map(
                        writer_process,
                        *zip(*proc_reads.items()),
                        itertools.repeat(log_queue),
                    )
                )

    log_queue.put("STOP")
    log_thread.join()

    logger.info("done!")
