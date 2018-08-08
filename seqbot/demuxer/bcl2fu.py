#!/usr/bin/env python

import glob
import gzip
import io
import itertools
import os
import struct

from collections import defaultdict, namedtuple, Counter

import numpy as np


CBCL = namedtuple('CBCL', ('version', 'header_size', 'bits_per_basecall',
                           'bits_per_qscore', 'number_of_bins', 'bins',
                           'number_of_tile_records', 'tile_offsets',
                           'non_PF_clusters_excluded'))

get_cycle = lambda cfn: int(os.path.basename(os.path.dirname(cfn))[1:-2])
get_part = lambda cfn: int(os.path.basename(cfn)[2])
get_tile = lambda cfn: int(os.path.basename(cfn)[4:8])



def read_cbcl_data(cbcl_files):
    cbcl_file_data = dict()

    for fn in cbcl_files:
        with open(fn, 'rb') as f:
            version, header_size, bits_per_basecall, bits_per_qscore, number_of_bins = struct.unpack(
                '<HIBBI', f.read(12))
            bins = np.fromfile(
                    f, dtype=np.uint32, count=2*number_of_bins
            ).reshape((number_of_bins, 2))

            number_of_tile_records = struct.unpack('<I', f.read(4))[0]
            tile_offsets = np.fromfile(
                f, dtype=np.uint32, count=4*number_of_tile_records
            ).reshape((number_of_tile_records, 4))

            non_PF_clusters_excluded = bool(struct.unpack('B', f.read(1))[0])

            cbcl_file_data[fn] = CBCL(version,
                                      header_size,
                                      bits_per_basecall,
                                      bits_per_qscore,
                                      number_of_bins,
                                      bins,
                                      number_of_tile_records,
                                      tile_offsets,
                                      non_PF_clusters_excluded)

    return cbcl_file_data


def read_cbcl_filters(filter_files):
    cbcl_filters = dict()

    for fn in filter_files:
        with open(fn, 'rb') as f:
            zv, filter_version, n_clusters = struct.unpack('<III', f.read(12))
            pf = np.fromfile(f, dtype=np.uint8, count=n_clusters)

            cbcl_filters[get_tile(fn)] = (pf & 0b1).astype(bool)

    return cbcl_filters


def cbcl_globber(bcl_path):
    cbcl_file_lists = dict()
    cbcl_filter_lists = dict()

    for lane in (1, 2, 3, 4):
        for part in itertools.count(1):
            cbcl_files = glob.glob(
                    os.path.join(bcl_path, 'Data', 'Intensities', 'BaseCalls',
                                 'L00{}'.format(lane), 'C*.1',
                                 'L00{}_{}.cbcl'.format(lane, part))
            )
            if cbcl_files:
                cbcl_files.sort(key=get_cycle)
                cbcl_file_lists[lane, part] = cbcl_files
            else:
                break

        cbcl_filter_list = glob.glob(
                os.path.join(bcl_path, 'Data', 'Intensities', 'BaseCalls',
                             'L00{}'.format(lane),
                             's_{}_*.filter'.format(lane))
        )

        cbcl_filter_list.sort(key=get_tile)
        cbcl_filter_lists[lane] = cbcl_filter_list

    return cbcl_file_lists, cbcl_filter_lists


def get_cbcl_data(cbcl_file_lists):
    cbcl_data = read_cbcl_data(cbcl_file_lists)

    number_of_tiles = {cbcl_data[fn].number_of_tile_records
                       for fn in cbcl_file_lists}

    assert len(number_of_tiles) == 1

    return cbcl_data, number_of_tiles.pop()


def get_byte_lists(cbcl_files, cbcl_data, cbcl_filter_data, tile_i):
    for fn in cbcl_files:
        ci = cbcl_data[fn]
        cf = cbcl_filter_data[ci.tile_offsets[tile_i, 0]]

        with open(fn, 'rb') as f:
            f.seek(ci.header_size + ci.tile_offsets[:tile_i, 3].sum(dtype=int))

            tile_data = io.BytesIO(f.read(ci.tile_offsets[tile_i, 3]))
            try:
                g = gzip.GzipFile(fileobj=tile_data, mode='r').read()
                byte_array = np.frombuffer(g, dtype=np.uint8,
                                           count=ci.tile_offsets[tile_i, 2])
            except OSError:
                yield None
                continue

            if ci.non_PF_clusters_excluded and cf.sum() % 2:
                yield np.hstack(((byte_array & 0b11)[:-1],
                                 (byte_array >> 4 & 0b11)))
            elif ci.non_PF_clusters_excluded:
                yield np.hstack(((byte_array & 0b11),
                                 (byte_array >> 4 & 0b11)))
            else:
                yield np.hstack(((byte_array & 0b11)[cf[::2]],
                                 (byte_array >> 4 & 0b11)[cf[1::2]]))


def extract_reads(cbcl_files, cbcl_filter_files, i, nproc):
    cbcl_filter_data = read_cbcl_filters(cbcl_filter_files)

    cbcl_data,n_tiles = get_cbcl_data(cbcl_files)

    for ii in range(i, n_tiles, nproc):
        ba_generator = enumerate(
            get_byte_lists(cbcl_files, cbcl_data, cbcl_filter_data, ii)
        )
        j, byte_array = next(ba_generator)

        byte_matrix = 4 * np.ones((byte_array.shape[0], len(cbcl_files)),
                                  dtype=np.uint8)
        byte_matrix[:, j] = byte_array

        for j, byte_array in ba_generator:
            if byte_array is not None:
                byte_matrix[:, j] = byte_array

        yield from (''.join('ACGTN'[b] for b in byte_matrix[k, :])
                    for k in range(byte_matrix.shape[0]))
