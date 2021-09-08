import re
import sys
import shlex
from typing import List
import subprocess as sp
from pathlib import Path


import numpy as np
import pandas as pd

from .tools import Plink, plink_exclude_across_bfiles
from .peakit import peakit
from .interactive_manh import output_peakplot

def read_assoc(filepath, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, logger, chunksize = 10000) -> pd.DataFrame:
    """
    Lazily load and filter the association file.
    
    Example
    -------
    >>> # GCTA file as input
    >>> assoc = read_assoc('/path/to/assoc.mlma.gz', 'Chr', 'bp', 'p', 'Freq', 'SNP', 'A1', 'A2')
    """
    assoc = pd.read_csv(filepath, sep = '\t',
                            chunksize = chunksize,
                            dtype = {
                                chr_col: np.int64,
                                pos_col: np.int64,
                                maf_col: np.float64,
                                rs_col: str,
                                a1_col: str,
                                a2_col: str
                            })
    for chunk in assoc:
        chunk[pval_col] = pd.to_numeric(chunk[pval_col], errors = 'coerce')
        nan_list = list(pd.isna(chunk[pval_col]))
        if nan_list.count(True):
            logger.debug(f"Removing {nan_list.count(True)} rows with invalid p-value")
        chunk = chunk.dropna().reset_index(drop = True)

        a1_check = chunk[a1_col].str.contains('[^ATGC]')
        if any(a1_check):
            invalid_strings = a1_check.to_list()
            logger.debug(f"Removing {invalid_strings.count(True)} rows with invalid a1 string value")
            logger.debug(chunk.loc[a1_check, [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]])
            chunk = chunk[~a1_check].reset_index(drop = True)

        a2_check = chunk[a2_col].str.contains('[^ATGC]')
        if any(a2_check):
            invalid_strings = a2_check.to_list()
            logger.debug(f"Removing {invalid_strings.count(True)} rows with invalid a2 string value")
            logger.debug(chunk.loc[a2_check, [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]])
            chunk = chunk[~a2_check].reset_index(drop = True)

        yield chunk


def get_signals(assoc, signif, chr_col, pos_col, pval_col) -> pd.DataFrame:
    """
    Accepts `read_assoc` function's iterable and outputs a DataFrame
    only including variants with p-value lower than the significance threshold.
    Outputs an empty DataFrame if no variants exceeds significance threshold.
    """
    concat_list = list()
    for chunk in assoc:
        chunk = chunk[signif>chunk[pval_col]]
        if chunk.shape[0]>0:
            concat_list.append(chunk)
    if len(concat_list)==0:
        return pd.DataFrame()
    signals = pd.concat(concat_list).reset_index(drop = True)
    signals.sort_values(by = [chr_col, pos_col], inplace = True)
    return signals


def run_locuszoom(build, peakdata_file, refsnp, rs_col, pval_col, db_file, prefix, ld, start, end, chrom):
    process = sp.run(shlex.split(f'''
        locuszoom \
            --build {build} \
            --metal {peakdata_file} \
            --refsnp "{refsnp}" \
            --markercol "{rs_col}" \
            --pvalcol "{pval_col}" \
            --db {db_file} \
            --prefix {prefix} \
            --plotonly showAnnot=T showRefsnpAnnot=T annotPch='"'21,24,24,25,22,22,8,7'"' rfrows=20 geneFontSize=.4 \
            --ld {ld} \
            --start={start} \
            --end={end} \
            --chr={chrom} showRecomb=T
        '''), stdout = sp.PIPE, stderr = sp.PIPE)
    return process


def _create_non_rs_to_pos_id(data: pd.DataFrame, chr_col, rs_col, pos_col) -> pd.Series:
    new_column = list()
    for chrom, _id, pos in data[[chr_col, rs_col, pos_col]].itertuples(index = False):
        if not _id.startswith('rs') and ':' not in _id:
            chrom = chrom.strip('chr')
            new_column.append(f'chr{chrom}:{pos}')
        else:
            new_column.append(_id)
    return pd.Series(new_column, dtype = str)


def _add_chr_to_id(column: pd.Series) -> pd.Series:
    new_column = list()
    for ele in column:
        if not ele.startswith('chr'):
            new_column.append(f'chr{ele}')
        else:
            new_column.append(ele)
    return pd.Series(new_column, dtype = str)


def _remove_build_suffix(column: pd.Series) -> pd.Series:
    """
    This is more of an ITG specific problem solving function.
    """
    return [re.sub('\[b3[78]\]', '', _id) for _id in column]

# TODO: Modulerise this function. Split it up to smaller functions
def process_peak(assocfile: str,
                  chr_col: str,
                  pos_col: str,
                  pval_col: str,
                  maf_col: str,
                  rs_col: str,
                  a1_col: str,
                  a2_col: str,
                  chrom: int,
                  start: int,
                  end: int,
                  current: int,
                  outdir: Path,
                  refflat: Path,
                  recomb: Path,
                  bfiles_list: List[str],
                  plink: Plink,
                  build: int,
                  ext_flank_kb: int,
                  logger):
    
    assoc = read_assoc(assocfile, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, logger)
    logger.info('Looking for signals..')
    concat_list = list()
    for chunk in assoc:
        filtered_chunk = chunk.loc[(chunk[chr_col]==chrom) & (start < chunk[pos_col]) & (chunk[pos_col] < end)]
        if filtered_chunk.shape[0] > 0:
            concat_list.append(filtered_chunk)
    logger.info('Generating peakdata')
    peakdata = pd.concat(concat_list).sort_values([chr_col, pos_col]).reset_index(drop = True)
    peakdata[rs_col] = _create_non_rs_to_pos_id(peakdata, chr_col, rs_col, pos_col)
    # '1:100' -> 'chr1:100'
    peakdata[rs_col] = _add_chr_to_id(peakdata[rs_col])
    # 'chr1:100[b38]' -> 'chr1:100'
    peakdata[rs_col] = _remove_build_suffix(peakdata[rs_col])
    peakdata_chrpos = peakdata[[rs_col, chr_col, pos_col]]
    peakdata_chrpos.columns = ['snp', 'chr', 'pos']

    peakdata_chrpos_path = outdir.joinpath(f'{chrom}.{start}.{end}.peakdata.chrpos')
    db_file = outdir.joinpath(f'{chrom}.{start}.db')
    logger.debug(f'Generating peakdata.chrpos to {peakdata_chrpos_path}')
    peakdata_chrpos.to_csv(peakdata_chrpos_path, sep = '\t', index = False)

    logger.info('Running dbmeister.py')
    logger.debug(sp.check_output(shlex.split(f"dbmeister.py --db {db_file} --snp_pos {peakdata_chrpos_path}")))
    logger.debug(sp.check_output(shlex.split(f"dbmeister.py --db {db_file} --refflat {refflat}")))
    logger.debug(sp.check_output(shlex.split(f"dbmeister.py --db {db_file} --recomb_rate {recomb}")))

    if start < 1:
        sensible_start = 1
        logger.warning(f'Negative start position changed to {sensible_start} : {chrom} {start} (1)')
    else:
        sensible_start = start

    logger.info(f'Extract {chrom}:{start}-{end} genotypes')
    mergelist = list()
    for count, bfile in enumerate(bfiles_list):
        out = outdir.joinpath(f'peak.{chrom}.{start}.{end}.{count}')
        logger.debug(f"plink.extract_genotypes('{bfile}', {chrom}, {start}, {end}, '{out}')")
        ps = plink.extract_genotypes(bfile, chrom, start, end, out)
        logger.debug(ps.stdout.decode())
        logger.debug(ps.stderr.decode())
        if ps.returncode==12:
            # TODO: Add some kind of test to ensure this part
            # A cohort can have no variants within a peak region
            # For example, peak was detected in chr1:200-300 driven by a variant
            # in cohortA, but cohortB has no variants within that region.
            logger.warning('Error Code 12 for plink.extract_genotypes')
            continue
        ## Modify BIM file 
        bimfile = f'{out}.bim'
        bim = pd.read_csv(bimfile, sep = '\t', header = None, names = ['chrom', 'id', '_', 'pos', 'a1', 'a2'])
        # 'non_rs_id' -> 'chr1:100'
        bim['id'] = _create_non_rs_to_pos_id(bim, 'chrom', 'id', 'pos')
        # '1:100' -> 'chr1:100'
        bim['id'] = _add_chr_to_id(bim['id'])
        # 'chr1:100[b38]' -> 'chr1:100'
        bim['id'] = _remove_build_suffix(bim['id'])
        bim.to_csv(bimfile, sep = '\t', header = False, index = False)
        mergelist.append(str(out))

    ## Merge plink binaries if multiple present
    assert len(mergelist) >= 1, f'mergelist length is {len(mergelist)}'
    logger.debug(f"{mergelist}")
    logger.info(f"Merging {mergelist}")
    mergelist_file = str(outdir.joinpath(f'{chrom}.{start}.{end}.mergelist'))
    out_merge = str(outdir.joinpath(f'{chrom}.{start}.{end}.peak'))
    
    logger.debug(f"plink.merge_region({mergelist_file}, '{mergelist}', {chrom}, {start}, {end}, '{out_merge}')")
    ps = plink.merge_region(mergelist_file, mergelist, chrom, start, end, out_merge)
    logger.debug(ps.stdout.decode())
    logger.debug(ps.stderr.decode())


    if ps.returncode == 3:
        logger.info('Multiallelic variants detected. Excluding')
        # Exclude variants which failed merge
        missnp_file = f'{out_merge}-merge.missnp'
        missnp_list = pd.read_csv(missnp_file, header = None)[0].to_list()
        peakdata = peakdata[~peakdata[rs_col].isin(missnp_list)].reset_index(drop = True)
        new_mergelist = plink_exclude_across_bfiles(plink, mergelist, missnp_file)
        
        # Make --merge-list file
        with open(mergelist_file, 'w') as f:
            for bfile in new_mergelist:
                f.write(f'{bfile}\n')
        
        # Delete old files
        for bfile in mergelist:
            Path(f'{bfile}.bed').unlink()
            Path(f'{bfile}.bim').unlink()
            Path(f'{bfile}.fam').unlink()
        logger.debug("plink.merge('{mergelist_file}', '{out_merge}')")
        ps = plink.merge(mergelist_file, out_merge)
        logger.debug(ps.stdout.decode())
        logger.debug(ps.stderr.decode())
        

    index_of_var_with_lowest_pval = peakdata[pval_col].idxmin()
    ref_snp_id = peakdata.loc[index_of_var_with_lowest_pval, rs_col]

    if ref_snp_id.startswith('rs'):
        refsnp = ref_snp_id # TODO: THIS CURRENTLY IS NEVER TRUE
    else:
        refsnp = ref_snp_id
        # chrom = str(peakdata.loc[index_of_var_with_lowest_pval, chr_col]).strip('chr')
        # pos = peakdata.loc[index_of_var_with_lowest_pval, pos_col]
        # refsnp = f'chr{chrom}:{pos}'

    logger.info(f"\n\nIn region {chrom} {start} {end}, top SNP is {refsnp}\n\n")

    logger.debug(f'plink.ld("{out_merge}", "{refsnp}", "{ext_flank_kb}", "{out_merge}")')
    ps = plink.ld(out_merge, refsnp, ext_flank_kb, out_merge)
    logger.debug(ps.stdout.decode())
    logger.debug(ps.stderr.decode())
    ld_data = pd.read_csv(f'{out_merge}.ld', delim_whitespace=True)
    ld_data = ld_data[['SNP_A', 'SNP_B', 'R2', 'R2']]
    ld_data.columns = ['snp1', 'snp2', 'dprime', 'rsquare']
    ld_file = f'{out_merge}.ld'
    ld_data.to_csv(ld_file, sep = ' ', header = True, index = False)

    subset_ld_data = ld_data[['snp2', 'dprime']].sort_values('snp2').reset_index(drop = True)
    subset_ld_data.columns = [rs_col, 'ld']

    peakdata_file = outdir.joinpath(f'{chrom}.{start}.{end}.peakdata.header')
    peakdata.to_csv(peakdata_file, sep = '\t', header = True, index = False)
    logger.debug('Merging LD info to main dataframe')
    joined_peakdata_ld = peakdata.merge(subset_ld_data, on = rs_col, how = 'left')
    joined_peakdata_ld['ld'].fillna(-0.01, inplace = True)
    joined_peakdata_ld['ld'] = joined_peakdata_ld['ld'].astype(float)
    logger.debug(joined_peakdata_ld['ld'])
    joined_peakdata_ld_file = outdir.joinpath(f'{chrom}.{start}.{end}.500kb')
    joined_peakdata_ld.to_csv(joined_peakdata_ld_file, sep = ',', header = True, index = False)

    logger.debug(f'run_locuszoom("{build}", "{peakdata_file}", "{refsnp}", "{rs_col}", "{pval_col}", "{db_file}", "{joined_peakdata_ld_file}", "{ld_file}", "{sensible_start}", "{end}", "{chrom}")')
    ps = run_locuszoom(build, peakdata_file, refsnp, rs_col, pval_col, db_file, joined_peakdata_ld_file, ld_file, sensible_start, end, chrom)
    logger.debug(ps.stdout.decode())
    logger.debug(ps.stderr.decode())
    if build == 38:
        b = 'b38'
    elif build == 37:
        b = 'b37'
    
    logger.debug("Running make_peakplot")
    input_file = str(joined_peakdata_ld_file)
    title = f"chr{chrom}:{start}-{end}"
    output_peakplot(
        infile = input_file,
        outfile = input_file, # Will output .html and .csv with this prefix
        title = title,
        pvalcol = pval_col,
        pscol = pos_col,
        mafcol = maf_col,
        chrcol = chr_col,
        a1col = a1_col,
        a2col = a2_col,
        build = b,
        logger = logger)

    logger.info(f"Done with peak {chrom} {start} {end}.")
    logger.info("Cleaning plink binary files")
    to_delete = list(outdir.glob(f'peak.{chrom}.{start}.{end}.*.*'))
    for file in to_delete:
        file.unlink()


def _make_done(outdir: Path):
    with open(outdir.joinpath('done'), 'w') as f:
        f.write('done')
    


def main(signif, assocfile, chr_col, pos_col, rs_col, pval_col, a1_col, a2_col, maf_col, bfiles, flank_bp, refflat, recomb, build, outdir, logger, memory = 30000):
    # ext_flank_bp = flank_bp + 100_000
    flank_kb = flank_bp // 1000
    ext_flank_kb = flank_kb + 100

    bfiles_list = bfiles.split(',')
    plink = Plink(memory)
    assoc = read_assoc(assocfile, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, logger)
    logger.debug('Getting signals')
    signals = get_signals(assoc, signif, chr_col, pos_col, pval_col)
    if signals.empty:
        logger.info("No peaks found. Exiting.")
        _make_done(outdir)
        sys.exit(0)
    logger.debug('Making peak_collections')
    peak_collections = peakit(signals, pval_col, chr_col, pos_col, flank_bp)
    peak_collections.merge()
    peaked = peak_collections.data
    logger.debug(peaked)
    peaked_file = outdir.joinpath('peaked')
    peaked.to_csv(peaked_file, sep = '\t', header = True, index = False)

    total_peak_count = peaked.shape[0]
    for current, (chrom, start, end) in peaked.iterrows():
        logger.info(f"Treating peak {chrom} {start} {end} (peak {current+1} / {total_peak_count} )")
        process_peak(assocfile,
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
                  current,
                  outdir,
                  refflat,
                  recomb,
                  bfiles_list,
                  plink,
                  build,
                  ext_flank_kb,
                  logger)
    _make_done(outdir)
    logger.info('Finished')
