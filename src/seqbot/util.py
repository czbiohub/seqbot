#!/usr/bin/env python

import csv
import gzip
import io
import itertools
import logging
import pathlib

from typing import Any, Dict, List, TextIO, Sequence

import concurrent.futures as cf

import boto3
import yaml
import xml.etree.ElementTree as et

import numpy as np


log = logging.getLogger(__name__)


default_config_file = pathlib.Path.home() / ".config/seqbot/config.yaml"


def get_config(config_file: pathlib.Path = default_config_file):
    log.debug(f"loading config_file {config_file}")
    with open(config_file) as f:
        config = yaml.load(f)

    return config


def read_samplesheet(samplesheet_io: TextIO):
    rows = list(csv.reader(samplesheet_io))

    # find the [Data] section to check format
    h_i = next(i for i, r in enumerate(rows) if r[0] == "[Data]")
    h_row = list(map(str.lower, rows[h_i + 1]))

    hdr = "\n".join(",".join(r) for r in rows[: h_i + 1])

    return hdr, h_row, rows[h_i + 2 :]


# hacky way to infer 10x version: check RunInfo.xml for the R1 length
def get_10x_version(seq_dir: pathlib.Path):
    with (seq_dir / "RunInfo.xml").open() as f:
        run_info = et.parse(f)

    read_elems = run_info.findall("./Run/Reads/Read[@NumCycles][@Number]")
    read_elems.sort(key=lambda el: int(el.get("Number")))

    r1_len = int(read_elems[0].get("NumCycles"))
    r2_len = int(read_elems[-1].get("NumCycles"))

    if r1_len == 26:
        return 2
    elif r1_len == 28 and len(read_elems) == 3:
        return 3
    elif r1_len == 28 and r2_len == 75 and len(read_elems) == 4:
        return "SR"  # for SpaceRanger although that might be too specific?
    elif r1_len == 150 and r2_len == 150:
        return "VDJ"
    else:
        return -1


def convert_index(
    row: List[str],
    index_i: int,
    index2_i: int,
    cr_indexes: Dict[str, Any],
    rev_comp: bool,
):
    if row[index_i] in cr_indexes:
        cr_index = cr_indexes[row[index_i]]
        if isinstance(cr_index, list):
            for cr_index in cr_indexes[row[index_i]]:
                yield [row[i] if i != index_i else cr_index for i in range(len(row))]
        else:
            new_row = [
                row[i] if i != index_i else cr_index["index1"] for i in range(len(row))
            ]

            index2 = cr_index["index2a"] if rev_comp else cr_index["index2b"]
            new_row.insert(index2_i, index2)

            yield new_row
    elif row[index_i].startswith("SI-"):
        raise ValueError(f"Unrecognized CellRanger index {row[index_i]}")
    else:
        yield row


def get_samplesheet(
    seq_dir: pathlib.Path, config: dict, download_path: pathlib.Path = None
):
    log.debug("creating S3 client")
    client = boto3.client("s3")

    log.info(f"downloading sample-sheet for {seq_dir.name}")

    fb = io.BytesIO()

    client.download_fileobj(
        Bucket=config["s3"]["seqbot_bucket"],
        Key=f"sample-sheets/{seq_dir.name}.csv",
        Fileobj=fb,
    )

    try:
        fb.getvalue().decode("ascii")
    except UnicodeDecodeError:
        log.error("illegal character in samplesheet")
        raise

    log.info(f"reading samplesheet for {seq_dir.name}")

    hdr, h_row, rows = read_samplesheet(
        io.StringIO(fb.getvalue().decode(), newline=None)
    )

    if "index" in h_row:
        index_i = h_row.index("index")
    else:
        raise ValueError("Samplesheet doesn't contain an index column")

    # if there's a lane column, we'll split lanes
    split_lanes = "lane" in h_row

    # hacky way to check for cellranger indexes:
    all_cellranger = all(r[index_i].startswith("SI-") for r in rows)
    any_cellranger = any(r[index_i].startswith("SI-") for r in rows)

    if any_cellranger and not all_cellranger:
        raise ValueError("Mix of CellRanger and other runs, can't demux this yet")

    if all_cellranger:
        with open(config["index"]["cellranger_indexes"]) as f:
            cr_indexes = {r[0]: r[1:] for r in csv.reader(f)}

        # new file for dual-indexed kits
        with open(config["index"]["cellranger_dual_indexes"]) as f:
            cr_indexes.update((r["index_name"], r) for r in csv.DictReader(f))

        cellranger = get_10x_version(seq_dir)
        if cellranger == "SR":
            index_i2 = index_i + 1
            h_row.insert(index_i2, "index2")
        else:
            index_i2 = None

        rev_comp = seq_dir.parent.name in config["seqs"]["rev-comp"]

        rows = [
            r
            for row in rows
            for r in convert_index(row, index_i, index_i2, cr_indexes, rev_comp)
        ]
    else:
        cellranger = 0

    if download_path is not None:
        with download_path.open("w") as out:
            print(hdr, file=out)
            print(",".join(h_row), file=out)
            print("\n".join(",".join(r) for r in rows), file=out)

    if "index2" in h_row:
        index2_i = h_row.index("index2")
        get_i = lambda r: r[index_i] + r[index2_i]
    else:
        get_i = lambda r: r[index_i]

    if split_lanes:
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


def hamming_conflict(indexes: Sequence[str], max_dist: int = 1):
    """Given a sequence of indexes, return True if any are within ``max_dist``
    of each other."""

    for d in range(max_dist + 1):
        h_sets = [hamming_set(index, d=d) for index in indexes]
        for hset1, hset2 in itertools.combinations(h_sets, 2):
            if hset1 & hset2:
                return True
    else:
        return False
