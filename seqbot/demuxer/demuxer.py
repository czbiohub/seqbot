#!/usr/bin/env python


# script to:
#   - scan sequencer folders
#   - check for completed runs
#   - demultiplex the run
#   - upload fastq files to S3 when completed


import csv
import glob
import io
import logging
import os
import pathlib
import smtplib
import subprocess
import sys
import time
import yaml

from email.mime.text import MIMEText
from logging.handlers import TimedRotatingFileHandler

import boto3

import utilities.log_util as ut_log
import utilities.s3_util as s3u


config_file = pathlib.Path('/home/seqbot/seqbot/config.yaml')

with open(config_file) as f:
    config = yaml.load(f)

# location where sequencer data gets written
SEQ_DIR = pathlib.Path(config['seqs']['base'])

# number of samples to batch when dealing with very large runs
sample_n = config['demux']['sample_n']

# cache files
local_samplesheets = pathlib.Path(config['cache']['samplesheet_dir'])
demux_cache = pathlib.Path(config['cache']['demuxed'])

def maybe_exit_process():
    # get all python pids
    pids = subprocess.check_output("pgrep python", shell=True).decode().split()

    cmds = 0
    for pid in pids:
        with open(pathlib.Path('/proc') / pid / 'cmdline') as f:
            # get the cmdline and match against this script
            line = f.read().split('\x00')[1]
            if line == sys.argv[0]:
                cmds += 1

    # if there are more than one, exit this one
    if cmds > 1:
        sys.exit(0)
    else:
        # otherwise, take a short nap before starting
        time.sleep(5)


def demux_mail(run_name:str):
    msg = MIMEText(
f'''Results are located in:
    s3://{config["s3"]["output_bucket"]}/{config["s3"]["fastq_prefix"]}/{run_name}

- seqbot
''')

    msg['Subject'] = f'[Seqbot] demux for {run_name} is complete!'
    msg['From'] = config['email']['username']
    msg['To'] = ','.join(config['email']['addresses_to_email'])

    with smtplib.SMTP('smtp.gmail.com', port=587) as smtp:
        smtp.ehlo()
        smtp.starttls()
        smtp.ehlo()
        smtp.login(config['email']['username'], config['email']['password'])

        smtp.send_message(msg)


def main(logger:logging.Logger, demux_set:set, samplesheets:set):
    logger.debug("Creating S3 client")
    client = boto3.client('s3')

    logger.info("Scanning {}...".format(SEQ_DIR))

    # for each sequencer, check for newly completed runs
    for seq in config['seqs']['dirs']:
        logger.info(seq)
        fns = list((SEQ_DIR / seq).glob(
                f'[0-9]*/{config["seqs"]["sentinels"][seq]}')
        )
        logger.debug(f'{len(fns)} ~complete runs in {SEQ_DIR / seq}')

        for fn in fns:
            seq_dir = fn.parent
            if seq_dir.name in demux_set:
                logger.debug(f'skipping {seq_dir.name}, already demuxed')
                continue

            if seq != 'NovaSeq-01' or (seq_dir / 'CopyComplete.txt').exists():
                if seq_dir.name in samplesheets:
                    logger.info(f'downloading sample-sheet for {seq_dir.name}')
                else:
                    logger.debug(f'skipping {seq_dir.name}, no sample-sheet')
                    continue

                fb = io.BytesIO()

                client.download_fileobj(
                    Bucket=config['s3']['samplesheet_bucket'],
                    Key=f'sample-sheets/{seq_dir.name}.csv',
                    Fileobj=fb
                )
                logger.info(f'reading samplesheet for {seq_dir.name}')
                rows = list(csv.reader(io.StringIO(fb.getvalue().decode())))

                # find the [Data] section to check format
                h_i = [i for i,r in enumerate(rows) if r[0] == '[Data]'][0]
                h_row = list(map(str.lower, rows[h_i + 1]))

                # if there's a lane column, we'll split lanes
                split_lanes = 'lane' in h_row

                if 'index' in h_row:
                    index_i = h_row.index('index')
                else:
                    logger.warning("Samplesheet doesn't contain an index column,"
                                " skipping!")
                    continue

                # hacky way to check for cellranger indexes:
                cellranger = rows[h_i + 2][index_i].startswith('SI-')

                # takes everything up to [Data]+1 line as header
                hdr = '\n'.join(','.join(r) for r in rows[:h_i + 2])
                rows = rows[h_i + 2:]
                batched = len(rows) > sample_n

                if cellranger and (split_lanes or batched):
                    logger.warning("Cellranger workflow won't use extra options")

                for i,j in enumerate(range(0, len(rows) + int(len(rows) % sample_n > 0), sample_n)):
                    with open(local_samplesheets / f'{seq_dir.name}_{i}.csv', 'w') as OUT:
                        print(hdr, file=OUT)
                        for r in rows[j:j + sample_n]:
                            print(','.join(r), file=OUT)

                for i in range(len(rows) + int(len(rows) % sample_n > 0) // sample_n):
                    logger.info(f'demuxing batch {i} of {seq_dir}')
                    demux_cmd = config['demux']['command_template'][:]
                    demux_cmd.extend(
                        ('-bcl_path',
                         f'local{seq_dir.as_uri()}',
                         '-output_path',
                         f's3://{config["s3"]["output_bucket"]}/{config["s3"]["fastq_prefix"]}/{seq_dir.name}',
                         '-samplesheet',
                         f'local{(local_samplesheets / seq_dir.name).as_uri()}_{i}.csv')
                    )

                    if not split_lanes:
                        demux_cmd.append('-no_lane_splitting')

                    if batched:
                        demux_cmd.extend(('-batch_runID', i+1))

                    if cellranger:
                        demux_cmd.append('-cellranger')

                    logger.debug(f"running command:\n\t{' '.join(demux_cmd)}")

                    subprocess.check_call(' '.join(demux_cmd), shell=True)

                logger.info('Sending notification email')
                demux_mail(seq_dir.name)
                demux_set.add(seq_dir.name)

    logger.info('scan complete')
    return demux_set


if __name__ == "__main__":
    # check for an existing process running
    maybe_exit_process()

    mainlogger = ut_log.get_trfh_logger(
        'demuxer',
        (config['logging']['info'], logging.INFO, 'W0', 10),
        (config['logging']['debug'], logging.DEBUG, 'midnight', 3)
    )

    if demux_cache.exists():
        mainlogger.debug('reading cache file for demuxed runs')
        with open(demux_cache) as f:
            demux_set = {line.strip() for line in f}
    else:
        mainlogger.debug('no cache file exists, querying S3...')
        demux_set = {
            fn.split('/', 2)[1] for fn in
            s3u.get_files(bucket=config['s3']['output_bucket'],
                          prefix=config['s3']['fastq_prefix'])
        }

    mainlogger.info(
        f'{len(demux_set)} folders in '
        f's3://{config["s3"]["output_bucket"]}/{config["s3"]["fastq_prefix"]}'
    )

    mainlogger.debug('Getting the list of sample-sheets')
    samplesheets = {
        os.path.splitext(os.path.basename(fn))[0] for fn in
        s3u.get_files(bucket=config['s3']['samplesheet_bucket'],
                      prefix='sample-sheets')
    }
    mainlogger.info(
        f'{len(samplesheets)} samplesheets in '
        f's3://{config["s3"]["samplesheet_bucket"]}/sample-sheets'
    )

    updated_demux_set = main(mainlogger, demux_set.copy(), samplesheets)

    mainlogger.info('demuxed {} new runs'.format(
            len(updated_demux_set) - len(demux_set))
    )
    with open(demux_cache, 'w') as OUT:
        print('\n'.join(sorted(updated_demux_set)), file=OUT)

    mainlogger.debug('wrote new cache file')
