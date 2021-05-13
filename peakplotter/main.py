#!/usr/bin/env python3
import sys

import click

from .errors import MissingExecutableError
from .utils import check_executable, DEPENDENT_EXECUTABLES


@click.command()
def cli(path):
    
    missing_executables = list()
    for exe in DEPENDENT_EXECUTABLES:
        if not check_executable(exe):
            missing_executables.append(exe)

    if missing_executables:
        raise MissingExecutableError(f"Executables missing: {', '.join(missing_executables)}")
        sys.exit(1)

if __name__ == '__main__':
    cli()