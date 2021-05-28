#!/usr/bin/env python3
import re

from setuptools import setup, find_packages


with open("peakplotter/__init__.py") as f:
    version = re.search(r"__version__ = \'(.*?)\'", f.read()).group(1)

setup(
    name = "peakplotter",
    version = version,
    install_requires = [
        'click',
        'numpy',
        'pandas',
        'bokeh',
        'requests',
        'pybedtools'
    ],
    include_package_data = True,
    entry_points = {
        'console_scripts': [
            'peakplotter = peakplotter.main:cli',
            'peakplotter-data-setup = peakplotter._data:setup_data',
        ],
    },
# Metadata
    author = "Arthur Gilly",
    author_email = "arthur.gilly@helmholtz-muenchen.de",
    packages = find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
)
