#!/usr/bin/env python3
"""
This script is for generating the centromere_b*.tsv files in peakplotter/data directory.
Regular users probably don't have any business using this file. 
"""


import shutil
import subprocess
from io import StringIO

import requests
import pandas as pd

def make_b38_centromere_data() -> pd.DataFrame:
    B38_CENTROMERE_URL = "https://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/38/Modeled_regions_for_GRCh38.tsv"

    s = requests.get(B38_CENTROMERE_URL).content
    cen_list = pd.read_csv(StringIO(s.decode('utf-8')), sep='\t', header=(0))
    cen_list.drop(['HET7'], inplace=True)
    cen_list.drop(cen_list.columns[3], axis=1, inplace=True)
    cen_list.columns = ['chrom','start','end']
    cen_list.drop(['CENX', 'CENY'], inplace = True)
    return cen_list

def make_b37_centromere_data() -> pd.DataFrame:
    B37_CENTROMERE_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz'
    cmd = f"wget -O- {B37_CENTROMERE_URL} | gunzip | grep -v -e chrX -e chrY | grep cen | mergeBed -i - | sed 's/chr//'"

    sp = subprocess.check_output(cmd, shell=True)

    cen_list = pd.read_csv(
                    StringIO(sp.decode('utf-8')),
                    sep = '\t',
                    header = None,
                    names = ['chrom', 'start', 'end']
    )

    cen_list = cen_list.sort_values('chrom').reset_index(drop = True)
    return cen_list


if __name__ == '__main__':
    missing_exe = list()
    for exe in ('wget', 'gunzip', 'grep', 'mergeBed', 'sed'):
        if shutil.which(exe) is None:
            missing_exe.append(exe)
    
    if missing_exe:
        raise RuntimeError(f"Following executables are missing in PATH: {', '.join(missing_exe)}")


    b38 = make_b38_centromere_data()
    b37 = make_b37_centromere_data()

    b38.to_csv('peakplotter/data/centromere_b38.tsv', sep = '\t', header = True, index = False)
    b37.to_csv('peakplotter/data/centromere_b37.tsv', sep = '\t', header = True, index = False)