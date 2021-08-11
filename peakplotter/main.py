#!/usr/bin/env python3
import sys
import shutil
import logging
from pathlib import Path

import click

from . import __version__
from ._data import get_data_path
from .utils import check_executable, DEPENDENT_EXECUTABLES
from .tools import Plink
from .errors import MissingExecutableError
from .plotpeaks import main, process_peak
from .logging import make_logger


@click.command()
@click.option('-a', '--assoc-file', type = click.Path(exists=True, dir_okay=False), required=True, help = 'Path to the association file. It can be gzipped, provided that it bears the .gz extension. Its first line must be a header, coherent with the name arguments below. It must be tab-separated, bgzipped and tabixed (tabix is available as part of bcftools)')
@click.option('-f', '--bfiles',type = click.STRING, required=True, help = 'Binary PLINK (.bed/.bim/.fam) file base name. This should contain the genotypes for at least all the variants in the assoc_file, but it can contain more. Please note that this is the base name, without the .bed/.bim/.fam extension.')
@click.option('-o', '--out', 'outdir', type = click.Path(file_okay=False, writable=True), required=True, help = 'Output directory to store all output files.')
@click.option('-chr', '--chr-col', type = click.STRING, required=True, help = 'Name of the column for chromosome names.')
@click.option('-ps', '--pos-col', type = click.STRING, required=True, help = 'Name of the column for chromosomal position.')
@click.option('-rs', '--rs-col', type = click.STRING, required=True, help = 'Name of the column for unique SNP ids (RS-id or chr:pos).')
@click.option('-p', '--pval-col', type = click.STRING, required=True, help = 'Name of the column for p-values.')
@click.option('-a1', '--a1-col', type = click.STRING, required=True, help = 'Name of the column for reference or major allele (used for predicting consequence).')
@click.option('-a2', '--a2-col', type = click.STRING, required=True, help = 'Name of the column for alternate or minor allele.')
@click.option('-maf', '--maf-col', type = click.STRING, required=True, help = 'Name of the column for non-reference or minor allele frequency.')
@click.option('-b', '--build', type = click.INT, default = 38, show_default=True, help = "Assembly build (37 or 38)")
@click.option('-s', '--signif', type=click.FLOAT, default=5e-8, help = 'The significance level above which to declare a variant significant. Scientific notation (such as 5e-8) is fine.')
@click.option('-bp', '--flank-bp', type = click.INT, default = 500_000, help = 'Flanking size in base pairs for drawing plots (defaults to 500kb, i.e. 1Mbp plots) around lead SNPs.')
@click.option('--debug', is_flag=True, flag_value = True, default = False, help = 'Set the log level from INFO to DEBUG.')
@click.option('--overwrite', is_flag=True, flag_value = True, default = False, help = 'Overwrite output directory if it already exists.')
@click.option('--version', is_flag=True, flag_value = True, default = False, help = 'Output version number of PeakPlotter')
def cli(assoc_file, bfiles, outdir, chr_col, pos_col, rs_col, pval_col, a1_col, a2_col, maf_col, build, signif, flank_bp, debug, overwrite, version):
    '''PeakPlotter
    '''
    if version is True:
        click.echo(__version__)
        sys.exit(0)


    ref_flat, recomb = get_data_path(build)
    if not ref_flat.exists() or not recomb.exists():
        click.echo('[ERROR] Download of some data is required. Please run peakplotter-data-setup in the commandline')
        sys.exit(1)

    # TODO: Change --bfile option type later to click.Path.
    # We do STRING for now, because the plotpeaks.sh script accepts comma-separated filepaths.
    missing_executables = list()
    for exe in DEPENDENT_EXECUTABLES:
        if not check_executable(exe):
            missing_executables.append(exe)

    if missing_executables:
        raise MissingExecutableError(f"Executables missing: {', '.join(missing_executables)}")
    
    
    outdir = Path(outdir)
    if outdir.exists() and overwrite is False:
        click.echo(f'[ERROR] Output directory ({str(outdir)}) already exists. Use --overwrite flag if you wish to overwrite.')
        sys.exit(1)
    if outdir.exists() and overwrite is True:
        click.echo(f'[INFO] Overwriting output directory ({str(outdir)}).')
        shutil.rmtree(outdir)
        outdir.mkdir()
    else:
        outdir.mkdir()

    # Save run configurations in the output log file
    configs = {
        'assoc_file': assoc_file,
        'bfiles': bfiles,
        'outdir': outdir,
        'chr_col': chr_col,
        'pos_col': pos_col,
        'rs_col': rs_col,
        'pval_col': pval_col,
        'a1_col': a1_col,
        'a2_col': a2_col,
        'maf_col': maf_col,
        'build': build,
        'signif': signif,
        'flank_bp': flank_bp
    }
    log_level = logging.DEBUG if debug else logging.INFO
    logger = make_logger(outdir.joinpath(f'{outdir.name}.log'), level = log_level)
    logger.info(f'PeakPlotter Version: {__version__}')
    arg_string = '\nArguments: \n'
    for k, v in configs.items():
        arg_string += f'  {k}: {v}\n'
    logger.info(arg_string)

    try:
        main(signif,
            assoc_file,
            chr_col,
            pos_col,
            rs_col,
            pval_col,
            a1_col,
            a2_col,
            maf_col,
            bfiles,
            flank_bp,
            ref_flat,
            recomb,
            build,
            outdir,
            logger,
            memory = 30000)
    except Exception:
        logger.critical('Unexpected error occurred:', exc_info = True)
        raise





@click.command()
@click.option('-a', '--assoc-file', type = click.Path(exists=True, dir_okay=False), required=True, help = 'Path to the association file. It can be gzipped, provided that it bears the .gz extension. Its first line must be a header, coherent with the name arguments below. It must be tab-separated, bgzipped and tabixed (tabix is available as part of bcftools)')
@click.option('-f', '--bfiles',type = click.STRING, required=True, help = 'Binary PLINK (.bed/.bim/.fam) file base name. This should contain the genotypes for at least all the variants in the assoc_file, but it can contain more. Please note that this is the base name, without the .bed/.bim/.fam extension.')
@click.option('-o', '--out', 'outdir', type = click.Path(file_okay=False, writable=True), required=True, help = 'Output directory to store all output files.')
@click.option('-chr', '--chr-col', type = click.STRING, required=True, help = 'Name of the column for chromosome names.')
@click.option('-ps', '--pos-col', type = click.STRING, required=True, help = 'Name of the column for chromosomal position.')
@click.option('-rs', '--rs-col', type = click.STRING, required=True, help = 'Name of the column for unique SNP ids (RS-id or chr:pos).')
@click.option('-p', '--pval-col', type = click.STRING, required=True, help = 'Name of the column for p-values.')
@click.option('-a1', '--a1-col', type = click.STRING, required=True, help = 'Name of the column for reference or major allele (used for predicting consequence).')
@click.option('-a2', '--a2-col', type = click.STRING, required=True, help = 'Name of the column for alternate or minor allele.')
@click.option('-maf', '--maf-col', type = click.STRING, required=True, help = 'Name of the column for non-reference or minor allele frequency.')
@click.option('-c', '--chrom', type = click.INT, required=True, help = 'Chromosome of the peak to plot')
@click.option('-s', '--start', type = click.INT, required=True, help = "Start of the peak to plot.")
@click.option('-e', '--end', type = click.INT, required=True, help = "End of the peak to plot.")
@click.option('-b', '--build', type = click.INT, default = 38, show_default=True, help = "Assembly build (37 or 38).")
@click.option('--debug', is_flag=True, flag_value = True, default = False, help = 'Set the log level from INFO to DEBUG.')
def cli_region(assoc_file, bfiles, outdir, chr_col, pos_col, rs_col, pval_col, a1_col, a2_col, maf_col, chrom, start, end, build, debug):

    ref_flat, recomb = get_data_path(build)
    if not ref_flat.exists() or not recomb.exists():
        click.echo('[ERROR] Download of some data is required. Please run peakplotter-data-setup in the commandline')
        sys.exit(1)

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
    
    outdir = Path(outdir)
    if not outdir.exists():
        outdir.mkdir()
    # Save run configurations in the output directory
    configs = {
        'run_mode': 'manual',
        'assoc_file': assoc_file,
        'bfiles': bfiles,
        'outdir': outdir,
        'chr_col': chr_col,
        'pos_col': pos_col,
        'rs_col': rs_col,
        'pval_col': pval_col,
        'a1_col': a1_col,
        'a2_col': a2_col,
        'maf_col': maf_col,
        'chrom': chrom,
        'start': start,
        'end': end,
        'build': build,
    }

    log_level = logging.DEBUG if debug else logging.INFO
    logger = make_logger(outdir.joinpath(f'{outdir.name}.{chrom}.{start}.{end}.log'), level = log_level)
    logger.info(f'PeakPlotter Version: {__version__}')
    arg_string = '\nArguments: \n'
    for k, v in configs.items():
        arg_string += f'  {k}: {v}\n'
    logger.info(arg_string)


    flank_bp = end - start
    flank_kb = flank_bp // 1000
    ext_flank_kb = flank_kb + 100

    bfiles_list = bfiles.split(',')
    plink = Plink(30_000)

    try:
        process_peak(assoc_file,
                  chr_col,
                  pos_col,
                  pval_col,
                  maf_col,
                  rs_col,
                  a1_col,
                  a2_col,
                  chrom,
                  start,
                  end,
                  1,
                  outdir,
                  ref_flat,
                  recomb,
                  bfiles_list,
                  plink,
                  build,
                  ext_flank_kb,
                  logger)
    except Exception:
        logger.critical('Unexpected error occurred:', exc_info = True)
        raise
    




if __name__ == '__main__':
    cli()
