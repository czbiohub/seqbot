#!/usr/bin/env python

import csv
import io
import itertools
import logging
import pathlib

import boto3
import yaml

import numpy as np


default_config_file = pathlib.Path.home() / ".config/seqbot/config.yaml"


def get_config(config_file: pathlib.Path = default_config_file):
    with open(config_file) as f:
        config = yaml.load(f)

    return config


def read_samplesheet(samplesheet_io): # should be typing.TextIO but this crashes??
    rows = list(csv.reader(samplesheet_io))

    # find the [Data] section to check format
    h_i = [i for i, r in enumerate(rows) if r[0] == "[Data]"][0]
    h_row = list(map(str.lower, rows[h_i + 1]))

    hdr = "\n".join(",".join(r) for r in rows[: h_i + 2])

    return hdr, h_row, rows[h_i + 2 :]


def get_samplesheet(seq_dir: pathlib.Path, config: dict, logger: logging.Logger):
    logger.debug("creating S3 client")
    client = boto3.client("s3")

    logger.info(f"downloading sample-sheet for {seq_dir.name}")

    fb = io.BytesIO()

    client.download_fileobj(
        Bucket=config["s3"]["seqbot_bucket"],
        Key=f"sample-sheets/{seq_dir.name}.csv",
        Fileobj=fb,
    )
    logger.info(f"reading samplesheet for {seq_dir.name}")

    hdr, h_row, rows = read_samplesheet(io.StringIO(fb.getvalue().decode(), newline=None))

    if "index" in h_row:
        index_i = h_row.index("index")
    else:
        logger.error("Samplesheet doesn't contain an index column, skipping!")
        raise ValueError("bad samplesheet")

    # if there's a lane column, we'll split lanes
    split_lanes = "lane" in h_row

    # hacky way to check for cellranger indexes:
    cellranger = rows[0][index_i].startswith("SI-")

    return hdr, rows, split_lanes, cellranger


def hamming_set(index: str, d: int = 1, include_N: bool = True):
    """Given an index of bases in {ACGTN}, generate all indexes within hamming
    distance d of the input

    :param index: string representing the index sequence
    :param d: maximum distance to allow
    :param include_N: include N when generating possible indexes
    :return: set of indexes within hamming distance d
    """

    base_d = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}

    new_base = [i * np.eye(len(index), dtype=np.uint8) for i in range(5 - include_N)]
    other_bases = 1 - np.eye(len(index), dtype=np.uint8)

    h_set = {tuple(base_d[c] for c in index)}

    for _ in range(d):
        for a in list(map(np.array, h_set)):
            h_set.update(t for i in new_base for t in map(tuple, a * other_bases + i))

    h_set = {"".join("ACGTN"[i] for i in h) for h in h_set}

    return h_set
