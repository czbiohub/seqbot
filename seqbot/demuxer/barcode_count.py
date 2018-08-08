#!/usr/bin/env python

import argparse
import gzip
import itertools
import logging
import os
import threading

from collections import defaultdict, Counter

import multiprocessing as mp

import seqbot.demuxer.bcl2fu as bcl2fu

import utilities.logging as ut_log


def get_parser():
    parser = argparse.ArgumentParser(
            prog='barcode_count.py',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--loglevel', type=int, default=logging.DEBUG)
    parser.add_argument('--n_threads', type=int, default=mp.cpu_count())

    parser.add_argument('--bcl_path', required=True)
    parser.add_argument('--output_dir', required=True)

    parser.add_argument('--index_cycle_start', required=True, type=int)
    parser.add_argument('--index_cycle_end', required=True, type=int)

    return parser


def read_count_processor(args):
    cbcl_files, cbcl_filter_files, i, nproc, out_file = args

    try:
        msg = 'starting pooljob with args: ({}..., {}..., {}, {})'.format(
            cbcl_files[0], cbcl_filter_files[0], i, nproc
        )
        log_queue.put((msg, logging.DEBUG))

        read_counter = Counter(read for read in bcl2fu.extract_reads(
            cbcl_files, cbcl_filter_files, i, nproc
        ))

        msg = 'pooljob done for args: ({}..., {}..., {}, {})'.format(
            cbcl_files[0], cbcl_filter_files[0], i, nproc
        )
        log_queue.put((msg, logging.DEBUG))

        log_queue.put(('writing to {}'.format(out_file), logging.INFO))
        with gzip.open(out_file, 'w') as OUT:
            for index in read_counter:
                OUT.write(
                    '{}\t{}\n'.format(index, read_counter[index]).encode()
                )

    except Exception as detail:
        log_queue.put(("encountered exception in process:\n{}".format(detail),
                       logging.INFO))


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    logger.setLevel(args.loglevel)

    cbcl_file_lists, cbcl_filter_lists = cbcl_globber(args.bcl_path)

    in_range = lambda cfn: (args.index_cycle_start
                            <= get_cycle(cfn)
                            < args.index_cycle_end)

    cbcl_file_lists = {
        (lane, part):tuple(cfn for cfn in cbcl_file_lists[lane, part]
                           if in_range(cfn))
        for lane,part in cbcl_file_lists
    }

    logger.info('{} CBCL files to read'.format(
            sum(map(len, cbcl_file_lists.values())))
    )

    lane_parts = sorted(cbcl_file_lists)

    global log_queue
    log_queue, log_thread = ut_log.get_thread_logger(logger)

    logger.debug('initializing pool of {} processes'.format(args.n_threads))

    pool = mp.Pool(args.n_threads)

    logger.info('reading {} files and aggregating counters'.format(
            sum(map(len, cbcl_file_lists.values()))
    ))

    output_file = os.path.join(args.output_dir, 'index_counts_{}.txt.gz')

    # warning: gratuitous use of itertools module ahead! it's gonna be great

    # lambda function to make this crazy itertools chain.
    # looks nuts, just repeats each element of s for [args.n_threads] times
    rep_n = lambda s: itertools.chain.from_iterable(
        map(itertools.repeat, s, itertools.repeat(args.n_threads))
    )

    # using imap_unordered to (maybe) keep memory usage low in the main thread
    try:
        logger.debug('starting count')
        for i,_ in enumerate(pool.imap_unordered(
            read_count_processor,
            zip(
                rep_n(cbcl_file_lists[lane, part] for lane,part in lane_parts),
                rep_n(cbcl_filter_lists[lane] for lane,part in lane_parts),
                itertools.cycle(range(args.n_threads)),
                itertools.repeat(args.n_threads),
                map(output_file.format, itertools.count())
            )
        )):
            if i % 100 == 0:
                logger.info(f'{i}')
    finally:
        pool.close()
        pool.join()

    log_queue.put('STOP')
    log_thread.join()

    logger.info('done!')


if __name__ == "__main__":
    mainlogger, log_file, file_handler = ut_log.get_logger('count_barcodes')

    try:
        main(mainlogger)
    except:
        mainlogger.info("An exception occurred", exc_info=True)
        raise
    finally:
        file_handler.close()
