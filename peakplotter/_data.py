"""
Module which contains code for downloading refFlat and recombination data files to run PeakPlotter.
"""

import os
import sys
import gzip
from io import BytesIO
from pathlib import Path

import click
import requests
import pandas as pd

from .data import REFFLAT_B37_PATH, REFFLAT_B38_PATH, RECOMB_B37_PATH, RECOMB_B38_PATH, DATA_DIR

_DOWNLOAD_URLS = {
    'HG19_REFFLAT_URL': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz',
    'HG19_RECOMB_URL': 'https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz',
    'HG38_REFFLAT_URL': 'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz',
    'HG38_RECOMB_URL': 'https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz',
}

def _download_refflat(build):
    if build==37:
        url = _DOWNLOAD_URLS['HG19_REFFLAT_URL']
    elif build==38:
        url = _DOWNLOAD_URLS['HG38_REFFLAT_URL']
    else:
        raise ValueError('Only 38 or 37 value allowed.')

    r = requests.get(url)
    bytes_data = gzip.decompress(r.content)
    column_names = ['geneName', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds']
    table = pd.read_table(BytesIO(bytes_data), sep = '\t', names = column_names)
    return table


def _download_recomb(build):
    if build==37:
        url = _DOWNLOAD_URLS['HG19_RECOMB_URL']
    elif build==38:
        url = _DOWNLOAD_URLS['HG38_RECOMB_URL']
    else:
        raise ValueError('Only 38 or 37 value allowed.')
    r = requests.get(url)
    bytes_data = gzip.decompress(r.content)
    table = pd.read_table(
                        BytesIO(bytes_data),
                        sep = ' ',
                        skiprows = 1,
                        names = ['chr', 'pos', 'recomb', 'cm_pos'])
    return table


def download_data(outdir):
    outdir_path = Path(outdir)
    if not outdir_path.exists() or not outdir_path.is_dir():
        raise NotADirectoryError(f"Either {outdir} doesn't exist or is not a directory.")
    
    writable = os.access(outdir_path, os.W_OK)
    if not writable:
        raise PermissionError(f"User does not have writable permission for {outdir}")
    
    to_csv_options = {'sep': '\t', 'index': False, 'doublequote': False}
    data = _download_refflat(37)
    data.to_csv(
        outdir_path.joinpath('refFlat_b37.tsv'),
        **to_csv_options
        )
    data = _download_refflat(38)
    data.to_csv(
        outdir_path.joinpath('refFlat_b38.tsv'),
        **to_csv_options
        )


    data = _download_recomb(37)
    data.to_csv(
        outdir_path.joinpath('recomb_rate_b37.tsv'),
        **to_csv_options
        )
    data = _download_recomb(38)
    data.to_csv(
        outdir_path.joinpath('recomb_rate_b38.tsv'),
        **to_csv_options
        )


@click.command()
@click.option('--force', is_flag=True, flag_value = True, default = False, help = 'Overwrite data directory even if it exists.')
def setup_data(force):
    """PeakPlotter data setup commandline interface
    This CLI is for downloading the necessary refFlat and recombination data to run PeakPlotter.
    """
    data_exists = [
        REFFLAT_B37_PATH.exists(),
        REFFLAT_B38_PATH.exists(),
        RECOMB_B37_PATH.exists(),
        RECOMB_B38_PATH.exists()
    ]
    
    if any(data_exists) and force is True:
        print(f'Overwriting data in {str(DATA_DIR)}')
        download_data(DATA_DIR)
    elif any(data_exists) and force is False:
        click.echo(f'At least one data in {str(DATA_DIR)} already exists.')
        click.echo('Use peakplotter-data-setup --force to overwrite data')
        sys.exit(0)
    else:
        click.echo(f'Downloading data to {str(DATA_DIR)}')
        download_data(DATA_DIR)


def get_data_path(build: int):
    if build==37:
        return REFFLAT_B37_PATH, RECOMB_B37_PATH
    if build==38:
        return REFFLAT_B38_PATH, RECOMB_B38_PATH
    else:
        raise ValueError('Only 38 or 37 value allowed.')