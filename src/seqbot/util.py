#!/usr/bin/env python

import csv
import gzip
import io
import itertools
import logging
import pathlib
import typing

import concurrent.futures as cf

import boto3
import yaml

import numpy as np


default_config_file = pathlib.Path.home() / ".config/seqbot/config.yaml"


def get_config(config_file: pathlib.Path = default_config_file):
    with open(config_file) as f:
        config = yaml.load(f)

    return config


def read_samplesheet(samplesheet_io: typing.TextIO):
    rows = list(csv.reader(samplesheet_io))

    # find the [Data] section to check format
    h_i = [i for i, r in enumerate(rows) if r[0] == "[Data]"][0]
    h_row = list(map(str.lower, rows[h_i + 1]))

    hdr = "\n".join(",".join(r) for r in rows[: h_i + 2])

    return hdr, h_row, rows[h_i + 2 :]


def get_samplesheet(
    seq_dir: pathlib.Path,
    config: dict,
    logger: logging.Logger,
    download_path: pathlib.Path = None,
):
    logger.debug("creating S3 client")
    client = boto3.client("s3")

    logger.info(f"downloading sample-sheet for {seq_dir.name}")

    fb = io.BytesIO()

    client.download_fileobj(
        Bucket=config["s3"]["seqbot_bucket"],
        Key=f"sample-sheets/{seq_dir.name}.csv",
        Fileobj=fb,
    )

    try:
        fb.getvalue().decode("ascii")
    except UnicodeDecodeError:
        logger.error("illegal character in samplesheet")
        raise

    logger.info(f"reading samplesheet for {seq_dir.name}")

    hdr, h_row, rows = read_samplesheet(
        io.StringIO(fb.getvalue().decode(), newline=None)
    )

    if download_path is not None:
        with download_path.open("w") as out:
            print(fb.getvalue().decode(), file=out)

    if "index" in h_row:
        index_i = h_row.index("index")
    else:
        logger.error("Samplesheet doesn't contain an index column, skipping!")
        raise ValueError("bad samplesheet")

    # if there's a lane column, we'll split lanes
    split_lanes = "lane" in h_row

    # hacky way to check for cellranger indexes:
    cellranger = any(r[index_i].startswith("SI-") for r in rows)

    if "index2" in h_row:
        index2_i = h_row.index("index2")
        get_i = lambda r: r[index_i] + r[index2_i]
    else:
        get_i = lambda r: r[index_i]

    if cellranger:
        # could still happen but we can't check it here
        index_overlap = False
    elif split_lanes:
        lane_i = h_row.index("lane")

        index_overlap = False

        for lane in {r[lane_i] for r in rows}:
            indexes = [get_i(r) for r in rows if r[lane_i] == lane]

            index_overlap = index_overlap or hamming_conflict(indexes, max_dist=1)
    else:
        indexes = [get_i(r) for r in rows]

        index_overlap = hamming_conflict(indexes, max_dist=1)

    return hdr, rows, split_lanes, cellranger, index_overlap


def hamming_set(index: str, d: int = 1, include_N: bool = True):
    """Given an index of bases in {ACGTN}, generate all indexes within hamming
    distance d of the input

    :param index: string representing the index sequence
    :param d: maximum distance to allow
    :param include_N: include N when generating possible indexes
    :return: set of indexes within hamming distance d
    """

    base_d = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}

    new_base = [i * np.eye(len(index), dtype=np.uint8) for i in range(4 + include_N)]
    other_bases = 1 - np.eye(len(index), dtype=np.uint8)

    h_set = {tuple(base_d[c] for c in index)}

    for _ in range(d):
        for a in list(map(np.array, h_set)):
            h_set.update(t for i in new_base for t in map(tuple, a * other_bases + i))

    h_set = {"".join("ACGTN"[i] for i in h) for h in h_set}

    return h_set


def hamming_conflict(indexes: typing.Sequence[str], max_dist: int = 1):
    """Given a sequence of indexes, return True if any are within ``max_dist``
    of each other."""

    for d in range(max_dist + 1):
        h_sets = [hamming_set(index, d=d) for index in indexes]
        for hset1, hset2 in itertools.combinations(h_sets, 2):
            if hset1 & hset2:
                return True
    else:
        return False


def write_file(seq_dir, fastq_name, fastqs):
    for fastq in fastqs:
        with gzip.open(seq_dir / fastq_name, "at") as out:
            with gzip.open(fastq, "rt") as f:
                for line in f:
                    out.write(line)


def merge_fastqs(seq_dir, n_processes, file_lists):
    with cf.ProcessPoolExecutor(max_workers=n_processes) as ex:
        list(ex.map(write_file, itertools.repeat(seq_dir), *zip(*file_lists.items())))
