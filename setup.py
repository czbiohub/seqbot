#!/usr/bin/env python

import glob
import os

import setuptools

version = "0.2"


setuptools.setup(
    name="seqbot",
    version=version,
    description="Scripts for sequencing automation",
    author="James Webber",
    author_email="james.webber@czbiohub.org",
    url="https://github.com/czbiohub/seqbot",
    packages=setuptools.find_packages(),
    install_requires=[
        "awscli >= 1.15.41",
        "awscli-cwlogs >= 1.4.4",
        "boto3 >= 1.7.72",
        "click >= 6.7",
        "PyYAML >= 3.12",
        "numpy >= 1.15.0",
        "utilities",
    ],
    long_description="See https://github.com/czbiohub/seqbot",
    license=open("LICENSE").readline().strip(),
    entry_points={
        "console_scripts": [
            "demuxer = seqbot.demuxer.demuxer:main",
            "nova_index = seqbot.demuxer.index_count:main",
            "nova_demux = seqbot.demuxer.write_fastq:main",
        ]
    },
)
