#!/usr/bin/env python3

import shlex
import subprocess

import click

from . import PLOTPEAKS_SCRIPT
from .errors import MissingExecutableError
from .utils import check_executable, _get_locuszoom_data_path, DEPENDENT_EXECUTABLES


@click.command()
@click.option('-a', '--assoc-file',type = click.Path(exists=True, dir_okay=False), help = 'Path to the association file. It can be gzipped, provided that it bears the .gz extension. Its first line must be a header, coherent with the name arguments below. It must be tab-separated, bgzipped and tabixed (tabix is available as part of bcftools)')
@click.option('-f', '--bfiles',type = click.STRING, help = 'Binary PLINK (.bed/.bim/.fam) file base name. This should contain the genotypes for at least all the variants in the assoc_file, but it can contain more. Please note that this is the base name, without the .bed/.bim/.fam extension.')
@click.option('-s', '--signif', type=click.FLOAT, help = 'The significance level above which to declare a variant significant. Scientific notation (such as 5e-8) is fine.')
@click.option('-chr', '--chr-col', type = click.STRING, help = 'Name of the column for chromosome names.')
@click.option('-ps', '--pos-col', type = click.STRING, help = 'Name of the column for chromosomal position.')
@click.option('-rs', '--rs-col', type = click.STRING, help = 'Name of the column for unique SNP ids (RS-id or chr:pos).')
@click.option('-p', '--pval-col', type = click.STRING, help = 'Name of the column for p-values.')
@click.option('-a1', '--a1-col', type = click.STRING, help = 'Name of the column for reference or major allele (used for predicting consequence).')
@click.option('-a2', '--a2-col', type = click.STRING, help = 'Name of the column for alternate or minor allele.')
@click.option('-maf', '--maf-col', type = click.STRING, help = 'Name of the column for non-reference or minor allele frequency.')
@click.option('-bp', '--flank-bp', type = click.INT, default = 500_000, help = 'Flanking size in base pairs for drawing plots (defaults to 500kb, i.e. 1Mbp plots) around lead SNPs.')
@click.option('--ref-flat', type = click.Path(exists=True, dir_okay=False), default = None, help = "Path to Locuszoom's refFlat file. By default, peakplotter finds it for you in the locuszoom files.")
@click.option('--recomb', type = click.Path(exists=True, dir_okay=False), default = None, help = "Path to Locuszoom's recomb_rate file. By default, peakplotter finds it for you in the locuszoom files.")
@click.option('--overwrite', is_flag=True, flag_value = True, default = False, help = 'Overwrite output directory if it already exists.')
def cli(assoc_file, bfiles, signif, chr_col, pos_col, rs_col, pval_col, a1_col, a2_col, maf_col, flank_bp, ref_flat, recomb, overwrite):
    '''PeakPlotter
    '''
    # TODO: Change --bfile option type later to click.Path.
    # We do STRING for now, because the plotpeaks.sh script accepts comma-separated filepaths.
    missing_executables = list()
    for exe in DEPENDENT_EXECUTABLES:
        if not check_executable(exe):
            missing_executables.append(exe)

    if missing_executables:
        raise MissingExecutableError(f"Executables missing: {', '.join(missing_executables)}")

    # TODO: Handle situation where only one file is given
    if (ref_flat is None) and (recomb is None):
        raise FileNotFoundError('Need to give ref_flat and recomb option')
        # ref_flat, recomb = _get_locuszoom_data_path() 

    command = f"{PLOTPEAKS_SCRIPT} {signif} {assoc_file} {chr_col} {pos_col} {rs_col} {pval_col} {a1_col} {a2_col} {maf_col} {bfiles} {flank_bp} {ref_flat} {recomb}"
    subprocess.run(shlex.split(command))

if __name__ == '__main__':
    cli()
