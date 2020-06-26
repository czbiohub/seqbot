#!/usr/bin/env python

import io
import glob
import os

import setuptools


def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ).read()


setuptools.setup(
    name="seqbot",
    version="0.3.1",
    license="MIT License",
    description="Scripts for sequencing automation",
    long_description=read("README.md"),
    author="James Webber",
    author_email="james.webber@czbiohub.org",
    url="https://github.com/czbiohub/seqbot",
    packages=setuptools.find_packages("src"),
    package_dir={"": "src"},
    py_modules=[
        os.path.splitext(os.path.basename(path))[0] for path in glob.glob("src/*.py")
    ],
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        "awscli >= 1.15.41",
        "awscli-cwlogs >= 1.4.4",
        "boto3 >= 1.7.72",
        "click >= 6.7",
        "PyYAML >= 3.12",
        "numpy >= 1.15.0",
        "czb-util",
    ],
    entry_points={
        "console_scripts": [
            "demuxer = seqbot.demuxer:main",
        ]
    },
)
