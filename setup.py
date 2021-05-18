#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name = "peakplotter",
    version = "0.0.2",
    install_requires = [
        'click',
        'numpy',
        'pandas',
        'bokeh',
        'requests'
    ],
    include_package_data = True,
    entry_points = {
        'console_scripts': [
            'peakplotter = peakplotter.main:cli',
            'peakplotter-data-setup = peakplotter.data:setup_data',
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
