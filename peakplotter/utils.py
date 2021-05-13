#!/usr/bin/env python3

import shutil

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




