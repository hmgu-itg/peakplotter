#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

from .errors import MissingExecutableError

DEPENDENT_EXECUTABLES = [
    'tabix', # HTSLib
    'bedtools', # BedTools
    'locuszoom', # Locuszoom
    'dbmeister.py', # Locuszoom
    'plink', # Plink
    'sponge' # Moreutils
]


def check_executable(executable: str) -> bool:
    '''
    Checks whether executable exists in PATH.

    Example
    -------
    >>> check_executable('cd')
    True
    >>> check_executable('non-existant-executable')
    False
    '''
    return bool(shutil.which(executable))



def _get_locuszoom_data_path() -> tuple[Path, Path]:
    locuszoom_path = shutil.which('locuszoom')
    if locuszoom_path is None:
        raise MissingExecutableError('locuszoom executable not in PATH')
    locuszoom_data = locuszoom_path.parent.parent.joinpath('data')
    recomb_file = locuszoom_data.joinpath('recomb_rate_b38.txt')
    ref_file = locuszoom_data.joinpath('refFlat_b38.txt')

    if not all([recomb_file.exists(), ref_file.exists()]):
        raise FileNotFoundError('Locuszoom data not found')

    return recomb_file, ref_file


