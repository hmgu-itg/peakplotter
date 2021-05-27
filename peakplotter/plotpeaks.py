import os
import shlex
from typing import List
import subprocess as sp
from pathlib import Path

import pandas as pd
import numpy as np

from .peakit import _peakit, bedtools_merge

def read_assoc(filepath, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, chunksize = 10000) -> pd.DataFrame:
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
                            pval_col: np.float64,
                            maf_col: np.float64,
                            rs_col: str,
                            a1_col: str,
                            a2_col: str
                        })

    return assoc


def get_signals(assoc, signif) -> pd.DataFrame:
    concat_list = list()
    for chunk in assoc:
        chunk = chunk[signif>chunk['p']]
        if chunk.shape[0]>0:
            concat_list.append(chunk)
    signals = pd.concat(concat_list).reset_index(drop = True)
    signals.sort_values(by = [chr_col, pos_col], inplace = True)
    return signals


class Plink:
    def __init__(self, memory = 30000):
        self.memory = memory
    
    @property
    def _cmd(self):
        return f'plink --memory {self.memory}'
    
    def __call__(self, command):
        return sp.run(shlex.split(f'{self._cmd} {command}'), capture_output=True)
    
    def dry_call(self, command):
        return f'{self._cmd} {command}'
    
    def extract_genotypes(self, bfile, chrom, start, end, out):
        return sp.run(shlex.split(f'{self._cmd} --bfile {bfile} --chr {chrom} --from-bp {start} --to-bp {end} --out {out} --make-bed'), capture_output=True)
    
    def merge(self, file: str, bfiles: List[str], chrom, start, end, out: str):
        with open(file, 'w') as f:
            for bfile in bfiles:
                f.write(f'{bfile}\n')
        process = sp.run(shlex.split(f'{self._cmd} --merge-list {file} --chr {chrom} --from-bp {start} --to-bp {end} --out {out} --make-bed'), capture_output=True)
        os.remove(file)
        return process
    
    def exclude(self, bfile, exclude, out):
        return sp.run(shlex.split(f'{self._cmd} --allow-no-sex --bfile {bfile} --exclude {exclude} --make-bed --out {out}'), capture_output=True)
    
    def ld(self, bfile, ld_snp, ext_flank_kb, out):
        return sp.run(shlex.split(f'{self._cmd} --allow-no-sex --bfile {bfile} --r2 --ld-snp {ld_snp}  --ld-window-kb {ext_flank_kb} --ld-window 999999 --ld-window-r2 0 --out {out}'), capture_output=True)


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
        '''), capture_output = True)
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
    

def process_peak(assocfile,
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
                  total_peak_count,
                  outdir,
                  refflat,
                  recomb,
                  bfiles_list,
                  plink,
                  build):
    print(f"Treating peak {chrom} {start} {end} (peak {current+1} / {total_peak_count} )")
    
    assoc = read_assoc(assocfile, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col)
    
    concat_list = list()
    for chunk in assoc:
        filtered_chunk = chunk.loc[(chunk[chr_col]==chrom) & (start < chunk[pos_col]) & (chunk[pos_col] < end)]
        if filtered_chunk.shape[0] > 0:
            concat_list.append(filtered_chunk)
    peakdata = pd.concat(concat_list).sort_values([chr_col, pos_col]).reset_index(drop = True)
    peakdata_chrpos = peakdata[[rs_col, chr_col, pos_col]].copy()
    # Add 'chr' to variant ID name
    # e.g. '1:100:A:G' -> 'chr1:100:A:G'
    rows_with_no_chr = ~peakdata_chrpos[rs_col].str.startswith('chr')
    peakdata_chrpos.loc[rows_with_no_chr, rs_col] = 'chr'+peakdata_chrpos.loc[rows_with_no_chr, rs_col]
    peakdata_chrpos.columns = ['snp', 'chr', 'pos']
    
    peakdata_chrpos_path = outdir.joinpath('peakdata.chrpos')
    db_file = outdir.joinpath(f'{chrom}.{start}.db')

    peakdata_chrpos.to_csv(peakdata_chrpos_path, sep = '\t', index = False)

    sp.check_output(shlex.split(f"dbmeister.py --db {db_file} --snp_pos {peakdata_chrpos_path}"))
    sp.check_output(shlex.split(f"dbmeister.py --db {db_file} --refflat {refflat}"))
    sp.check_output(shlex.split(f"dbmeister.py --db {db_file} --recomb_rate {recomb}"))
    
    if start < 1:
        sensible_start = 1
        print(f'\n\n\nWARNING\t Negative start position changed to {sensible_start} : {chrom} {start} (1)\n\n\n')
    else:
        sensible_start = start



    for count, bfile in enumerate(bfiles_list):
        out = outdir.joinpath(f'peak.{chrom}.{start}.{end}.{count}')
        print(plink.extract_genotypes(bfile, chrom, start, end, out))
        
        ## Modify BIM file 
        bimfile = f'{out}.bim'
        bim = pd.read_csv(bimfile, sep = '\t', header = None, names = ['chrom', 'id', '_', 'pos', 'a1', 'a2'])

        # 'non_rs_id' -> 'chr1:100'
        bim['id'] = _create_non_rs_to_pos_id(bim, 'chrom', 'id', 'pos')

        # '1:100' -> 'chr1:100'
        bim['id'] = _add_chr_to_id(bim['id'])
        bim.to_csv(bimfile, sep = '\t', header = False, index = False)
        
        

    ## Merge plink binaries if multiple present
    mergelist = [str(f).strip('.bed') for f in outdir.glob(f'peak.{chrom}.{start}.{end}.*.bed')]
    assert len(mergelist) >= 1, f'mergelist length is {len(mergelist)}'

    if len(mergelist)==1:
        bfile = Path(mergelist[0])
        out_merge = bfile.parent.joinpath('merged')
        Path(f'{bfile}.bed').rename(f'{out_merge}.bed')
        Path(f'{bfile}.bim').rename(f'{out_merge}.bim')
        Path(f'{bfile}.fam').rename(f'{out_merge}.fam')
    elif len(mergelist) > 1:
        mergelist_file = str(outdir.joinpath('tmp_mergelist'))
        out_merge = str(outdir.joinpath('merged'))
        
        process = plink.merge(mergelist_file, mergelist, chrom, start, end, out)
        if process.returncode == 3:
            # Exclude variants which failed merge
            missnp_file = f'{out}-merge.missnp'
            missnp_list = pd.read_csv(missnp_file, header = None)[0].to_list()
            peakdata = peakdata[~peakdata[rs_col].isin(missnp_list)].reset_index(drop = True)
            for bfile in mergelist:
                plink.exclude(bfile, missnp_file, f'{bfile}.tmp')
                Path(f'{bfile}.tmp.bed').rename(f'{bfile}.bed')
                Path(f'{bfile}.tmp.bim').rename(f'{bfile}.bim')
                Path(f'{bfile}.tmp.fam').rename(f'{bfile}.fam')
                
            process = plink.merge(mergelist_file, mergelist, chrom, start, end, out)
        
        
        
    peakdata[rs_col] = _create_non_rs_to_pos_id(peakdata, chr_col, rs_col, pos_col)
    peakdata[rs_col] = _add_chr_to_id(peakdata[rs_col])

    index_of_var_with_lowest_pval = peakdata[pval_col].idxmin()
    ref_snp_id = peakdata.loc[index_of_var_with_lowest_pval, rs_col]

    if ref_snp_id.startswith('rs'):
        refsnp = ref_snp_id
    else:
        chrom = str(peakdata.loc[index_of_var_with_lowest_pval, chr_col]).strip('chr')
        pos = peakdata.loc[index_of_var_with_lowest_pval, pos_col]
        refsnp = f'chr{chrom}:{pos}'

    print(f"\n\nIn region {chrom} {start} {end}, top SNP is {refsnp}\n\n")

    plink.ld(bfile, refsnp, ext_flank_kb, out_merge)
    ld_data = pd.read_csv(f'{out_merge}.ld', delim_whitespace=True)
    ld_data = ld_data[['SNP_A', 'SNP_B', 'R2', 'R2']]
    ld_data.columns = ['snp1', 'snp2', 'dprime', 'rsquare']
    ld_file = f'{out_merge}.ld'
    ld_data.to_csv(ld_file, sep = ' ', header = True, index = False)

    peakdata_file = outdir.joinpath('peakdata.header')
    peakdata.to_csv(peakdata_file, sep = '\t', header = True, index = False)


    run_locuszoom(build, peakdata_file, refsnp, rs_col, pval_col, db_file, f'{chrom}.{start}.{end}.500kb', ld_file, sensible_start, end, chrom)


def main(signif, assocfile, chr_col, pos_col, rs_col, pval_col, a1_col, a2_col, maf_col, bfiles, flank_bp, refflat, recomb, build, outdir, memory = 30000):
    ext_flank_bp = flank_bp + 100_000
    flank_kb = flank_bp // 1000
    ext_flank_kb = flank_kb + 100

    bfiles_list = bfiles.split(',')
    plink = Plink(memory)
    assoc = read_assoc(assocfile, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col)
    signals = get_signals(assoc, signif)
    peak_collections = _peakit(signals, pval_col, chr_col, pos_col)
    peaked = bedtools_merge(peak_collections.data)


    total_peak_count = peaked.shape[0]
    for current, (chrom, start, end) in peaked.iterrows():
        print(f"Treating peak {chrom} {start} {end} (peak {current+1} / {total_peak_count} )")
        
        assoc = read_assoc(assocfile, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col)
        
        concat_list = list()
        for chunk in assoc:
            filtered_chunk = chunk.loc[(chunk[chr_col]==chrom) & (start < chunk[pos_col]) & (chunk[pos_col] < end)]
            if filtered_chunk.shape[0] > 0:
                concat_list.append(filtered_chunk)
        peakdata = pd.concat(concat_list).sort_values([chr_col, pos_col]).reset_index(drop = True)
        peakdata_chrpos = peakdata[[rs_col, chr_col, pos_col]].copy()
        # Add 'chr' to variant ID name
        # e.g. '1:100:A:G' -> 'chr1:100:A:G'
        rows_with_no_chr = ~peakdata_chrpos[rs_col].str.startswith('chr')
        peakdata_chrpos.loc[rows_with_no_chr, rs_col] = 'chr'+peakdata_chrpos.loc[rows_with_no_chr, rs_col]
        peakdata_chrpos.columns = ['snp', 'chr', 'pos']
        
        peakdata_chrpos_path = outdir.joinpath('peakdata.chrpos')
        db_file = outdir.joinpath(f'{chrom}.{start}.db')

        peakdata_chrpos.to_csv(peakdata_chrpos_path, sep = '\t', index = False)

        sp.check_output(shlex.split(f"dbmeister.py --db {db_file} --snp_pos {peakdata_chrpos_path}"))
        sp.check_output(shlex.split(f"dbmeister.py --db {db_file} --refflat {refflat}"))
        sp.check_output(shlex.split(f"dbmeister.py --db {db_file} --recomb_rate {recomb}"))
        
        if start < 1:
            sensible_start = 1
            print(f'\n\n\nWARNING\t Negative start position changed to {sensible_start} : {chrom} {start} (1)\n\n\n')
        else:
            sensible_start = start



        for count, bfile in enumerate(bfiles_list):
            out = outdir.joinpath(f'peak.{chrom}.{start}.{end}.{count}')
            print(plink.extract_genotypes(bfile, chrom, start, end, out))
            
            ## Modify BIM file 
            bimfile = f'{out}.bim'
            bim = pd.read_csv(bimfile, sep = '\t', header = None, names = ['chrom', 'id', '_', 'pos', 'a1', 'a2'])

            # 'non_rs_id' -> 'chr1:100'
            bim['id'] = _create_non_rs_to_pos_id(bim, 'chrom', 'id', 'pos')

            # '1:100' -> 'chr1:100'
            bim['id'] = _add_chr_to_id(bim['id'])
            bim.to_csv(bimfile, sep = '\t', header = False, index = False)
            
            

        ## Merge plink binaries if multiple present
        mergelist = [str(f).strip('.bed') for f in outdir.glob(f'peak.{chrom}.{start}.{end}.*.bed')]
        assert len(mergelist) >= 1, f'mergelist length is {len(mergelist)}'

        if len(mergelist)==1:
            bfile = Path(mergelist[0])
            out_merge = bfile.parent.joinpath('merged')
            Path(f'{bfile}.bed').rename(f'{out_merge}.bed')
            Path(f'{bfile}.bim').rename(f'{out_merge}.bim')
            Path(f'{bfile}.fam').rename(f'{out_merge}.fam')
        elif len(mergelist) > 1:
            mergelist_file = str(outdir.joinpath('tmp_mergelist'))
            out_merge = str(outdir.joinpath('merged'))
            
            process = plink.merge(mergelist_file, mergelist, chrom, start, end, out)
            if process.returncode == 3:
                # Exclude variants which failed merge
                missnp_file = f'{out}-merge.missnp'
                missnp_list = pd.read_csv(missnp_file, header = None)[0].to_list()
                peakdata = peakdata[~peakdata[rs_col].isin(missnp_list)].reset_index(drop = True)
                for bfile in mergelist:
                    plink.exclude(bfile, missnp_file, f'{bfile}.tmp')
                    Path(f'{bfile}.tmp.bed').rename(f'{bfile}.bed')
                    Path(f'{bfile}.tmp.bim').rename(f'{bfile}.bim')
                    Path(f'{bfile}.tmp.fam').rename(f'{bfile}.fam')
                    
                process = plink.merge(mergelist_file, mergelist, chrom, start, end, out)
            
            
            
        peakdata[rs_col] = _create_non_rs_to_pos_id(peakdata, chr_col, rs_col, pos_col)
        peakdata[rs_col] = _add_chr_to_id(peakdata[rs_col])

        index_of_var_with_lowest_pval = peakdata[pval_col].idxmin()
        ref_snp_id = peakdata.loc[index_of_var_with_lowest_pval, rs_col]

        if ref_snp_id.startswith('rs'):
            refsnp = ref_snp_id
        else:
            chrom = str(peakdata.loc[index_of_var_with_lowest_pval, chr_col]).strip('chr')
            pos = peakdata.loc[index_of_var_with_lowest_pval, pos_col]
            refsnp = f'chr{chrom}:{pos}'

        print(f"\n\nIn region {chrom} {start} {end}, top SNP is {refsnp}\n\n")

        plink.ld(bfile, refsnp, ext_flank_kb, out_merge)
        ld_data = pd.read_csv(f'{out_merge}.ld', delim_whitespace=True)
        ld_data = ld_data[['SNP_A', 'SNP_B', 'R2', 'R2']]
        ld_data.columns = ['snp1', 'snp2', 'dprime', 'rsquare']
        ld_file = f'{out_merge}.ld'
        ld_data.to_csv(ld_file, sep = ' ', header = True, index = False)

        peakdata_file = outdir.joinpath('peakdata.header')
        peakdata.to_csv(peakdata_file, sep = '\t', header = True, index = False)


        run_locuszoom(build, peakdata_file, refsnp, rs_col, pval_col, db_file, f'{chrom}.{start}.{end}.500kb', ld_file, sensible_start, end, chrom)