#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name = "peakplotter",
    version = "0.0.1",
    install_requires = [
        'asr',
        'numpy',
        'pandas',
        'bokeh',
        'requests',

    ],
    include_package_data = False,
    
# Metadata
    author = "Arthur Gilly",
    author_email = "arthur.gilly@helmholtz-muenchen.de",
    packages = find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
)