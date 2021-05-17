"""
Module which contains code for downloading refFlat and recombination data files to run PeakPlotter.
"""

import os
import gzip
from io import BytesIO
from pathlib import Path

import click
import requests
import pandas as pd

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
    base = Path(__file__).absolute().parent
    data_dir = base.joinpath('data')
    if not data_dir.exists():
        print(f'Downloading data to {str(data_dir)}')
        data_dir.mkdir()
        download_data(data_dir)
    elif data_dir.exists() and force is True:
        print(f'Overwriting data to {str(data_dir)}')
        download_data(data_dir)
    else:
        print(f'{str(data_dir)} directory already exists.')


def get_data_path(build):
    base = Path(__file__).absolute().parent
    data_dir = base.joinpath('data')
    if build==37:
        return data_dir.joinpath('refFlat_b37.tsv'), data_dir.joinpath('recomb_b37.tsv')
    if build==38:
        return data_dir.joinpath('refFlat_b38.tsv'), data_dir.joinpath('recomb_b38.tsv')
    else:
        raise ValueError('Only 38 or 37 value allowed.')