import pandas as pd
import numpy as np

def read_assoc(filepath, signif, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, chunksize = 10000) -> pd.DataFrame:
    """
    Lazily load and filter the association file.
    
    Example
    -------
    >>> # GCTA file as input
    >>> signals = read_assoc('/path/to/assoc.mlma.gz', 5e-8, 'Chr', 'bp', 'p', 'Freq', 'SNP', 'A1', 'A2')
    """
    chunks = pd.read_csv(filepath, sep = '\t',
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

    concat_list = list()
    for chunk in chunks:
        chunk = chunk[signif>chunk['p']]
        if chunk.shape[0]>0:
            concat_list.append(chunk)
    data = pd.concat(concat_list).reset_index(drop = True)
    data.sort_values(by = [chr_col, pos_col], inplace = True)
    return data