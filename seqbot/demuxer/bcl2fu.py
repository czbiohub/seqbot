#!/usr/bin/env python

import glob
import gzip
import itertools
import pathlib
import struct

from collections import defaultdict, namedtuple, Counter

import numpy as np


CBCL = namedtuple('CBCL', ('version', 'header_size', 'bits_per_basecall',
                           'bits_per_qscore', 'number_of_bins', 'bins',
                           'number_of_tile_records', 'tile_offsets',
                           'non_PF_clusters_excluded'))

get_cycle = lambda cfn: int(cfn.parent.name[1:-2])
get_part = lambda cfn: int(cfn.name[2])
get_tile = lambda cfn: int(cfn.name[4:8])



def read_cbcl_data(cbcl_files):
    cbcl_data = dict()

    for fn in cbcl_files:
        with open(fn, 'rb') as f:
            version, header_size, bits_per_basecall, bits_per_qscore, number_of_bins = struct.unpack(
                '<HIBBI', f.read(12))
            bins = np.fromfile(
                f, dtype=np.uint32, count=2*number_of_bins
            ).reshape((number_of_bins, 2)).astype(np.uint8)

            number_of_tile_records = struct.unpack('<I', f.read(4))[0]
            tile_offsets = np.fromfile(
                f, dtype=np.uint32, count=4*number_of_tile_records
            ).reshape((number_of_tile_records, 4))

            non_PF_clusters_excluded = bool(struct.unpack('B', f.read(1))[0])

            cbcl_data[fn] = CBCL(version,
                                 header_size,
                                 bits_per_basecall,
                                 bits_per_qscore,
                                 number_of_bins,
                                 bins,
                                 number_of_tile_records,
                                 tile_offsets,
                                 non_PF_clusters_excluded)

    return cbcl_data


def read_cbcl_filters(filter_files):
    cbcl_filters = dict()

    for fn in filter_files:
        with open(fn, 'rb') as f:
            _, filter_version, n_clusters = struct.unpack('<III', f.read(12))
            pf = np.fromfile(f, dtype=np.uint8, count=n_clusters)

            cbcl_filters[get_tile(fn)] = (pf & 0b1).astype(bool)

    return cbcl_filters


def read_cbcl_locs(loc_file):
    with open(loc_file, 'rb') as f:
        _, _, n_clusters = struct.unpack('<III', f.read(12))
        cbcl_locs = (np.fromfile(
                f, dtype=np.float32, count=2 * n_clusters
        ) * 10 + 1000).reshape((n_clusters, 2)).round().astype(int)

    return cbcl_locs


def cbcl_id(bcl_path: pathlib.Path):
    _, instrument, run_id, flowcell_id = bcl_path.name.split('_')

    return f'@{instrument}:{int(run_id)}:{flowcell_id}'


def cbcl_globber(bcl_path: pathlib.Path):
    cbcl_file_lists = dict()
    cbcl_filter_lists = dict()

    for lane in (1, 2, 3, 4):
        for part in itertools.count(1):
            cbcl_files = sorted((
                bcl_path / 'Data' / 'Intensities' / 'BaseCalls' / f'L00{lane}'
            ).glob(f'C*.1/L00{lane}_{part}.cbcl'), key=get_cycle)

            if cbcl_files:
                cbcl_file_lists[lane, part] = cbcl_files
            else:
                break

        cbcl_filter_list = sorted((
            bcl_path / 'Data' / 'Intensities' / 'BaseCalls' / f'L00{lane}'
        ).glob(f's_{lane}_*.filter'), key=get_tile)

        cbcl_filter_lists[lane] = cbcl_filter_list

    loc_file = bcl_path / 'Data' / 'Intensities' / 's.locs'
    assert loc_file.exists()

    return cbcl_file_lists, cbcl_filter_lists, loc_file


def get_cbcl_data(cbcl_file_list):
    cbcl_data = read_cbcl_data(cbcl_file_list)

    # assertion: the tile numbers are the same across all of these files
    assert np.array_equal(
        np.vstack([ci.tile_offsets[:,0] for ci in cbcl_data.values()]).max(0),
        np.vstack([ci.tile_offsets[:,0] for ci in cbcl_data.values()]).min(0)
    )

    tiles = cbcl_data[cbcl_file_list[0]].tile_offsets[:,0]

    return cbcl_data, tiles


def get_byte_lists(cbcl_files, cbcl_data, cbcl_filter_data, tile_i):
    for fn in cbcl_files:
        ci = cbcl_data[fn]
        cf = cbcl_filter_data[ci.tile_offsets[tile_i, 0]]
        odd = -(cf.sum() % 2) or None

        with open(fn, 'rb') as f:
            f.seek(ci.header_size + ci.tile_offsets[:tile_i, 3].sum(dtype=int))

            try:
                byte_array = np.frombuffer(
                    gzip.decompress(f.read(ci.tile_offsets[tile_i, 3])),
                    dtype=np.uint8, count=ci.tile_offsets[tile_i, 2]
                )
            except OSError:
                yield (None, None)
                continue

            ub = np.unpackbits(byte_array[::-1, None], axis=1).reshape((-1, 2))
            pb = (np.packbits(ub, axis=1) >> 6).squeeze()[::-1]

            base_array = pb[::2]
            qscore_array = pb[1::2]

            if not ci.non_PF_clusters_excluded:
                yield base_array[cf], qscore_array[cf]
            else:
                yield base_array[:odd], qscore_array[:odd]


def extract_reads(cbcl_files, cbcl_filter_files, loc_file, i, nproc):
    cbcl_filter_data = read_cbcl_filters(cbcl_filter_files)
    cbcl_data, tiles = get_cbcl_data(cbcl_files)
    cbcl_locs = read_cbcl_locs(loc_file)

    # retrieve qscore binning
    code_set = set(tuple(ci.bins[:,1] + 33) for ci in cbcl_data.values())

    # assuming for now that bins are constant across a run
    assert len(code_set) == 1
    code_set = bytes(code_set.pop()).decode() + '#'

    for ii in range(i, tiles.shape[0], nproc):
        tile_no = tiles[ii]
        cf = cbcl_filter_data[tile_no]

        ba_generator = enumerate(
            get_byte_lists(cbcl_files, cbcl_data, cbcl_filter_data, ii)
        )
        j, (base_array, qscore_array) = next(ba_generator)

        base_matrix = 4 * np.ones((base_array.shape[0], len(cbcl_files)),
                                  dtype=np.uint8)
        base_matrix[:, j] = base_array
        qscore_matrix = 4 * np.ones_like(base_matrix)
        qscore_matrix[:, j] = qscore_array

        for j, (base_array, qscore_array) in ba_generator:
            if base_array is not None:
                base_matrix[:, j] = base_array
                qscore_matrix[:, j] = qscore_array

        base_matrix[qscore_matrix == 0] = 4
        qscore_matrix[qscore_matrix == 0] = 4

        yield from zip((''.join('ACGTN'[b] for b in base_matrix[k, :])
                        for k in range(base_matrix.shape[0])),
                       (''.join(code_set[b] for b in qscore_matrix[k, :])
                        for k in range(qscore_matrix.shape[0])),
                       (f'{tile_no}:{r[0]}:{r[1]}' for r in cbcl_locs[cf,:]))
