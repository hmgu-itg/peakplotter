import io 

import pandas as pd

from peakplotter import plotpeaks
from peakplotter.test_utils import get_test_logger

def test_read_assoc_with_METAL_output():
    """
    METAL outputs have 'nanenan' values in p-value columns. Exclude these.
    
    Note
    ----
    1. https://genome.sph.umich.edu/wiki/METAL_Documentation
    """
    chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col = 'c', 'pos', 'pval', 'maf', 'rs', 'a1', 'a2'
    example = pd.DataFrame([
            [1, 'rs001', 100, 'G', 'T', 0.01, 'nanenan'],
            [1, 'rs002', 200, 'G', 'T', 0.01, '0.01']
        ]
        , columns = [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]
    )

    example_buff = io.StringIO()
    example.to_csv(example_buff, sep = '\t', header = True, index = False)
    example_buff.seek(0)
    logger = get_test_logger()
    iterable = plotpeaks.read_assoc(example_buff, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, logger)
    data = next(iterable)

    expected = pd.DataFrame([
            [1, 'rs002', 200, 'G', 'T', 0.01, 0.01]
        ]
        , columns = [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]
    )
    
    pd.testing.assert_frame_equal(expected, data)



def test_read_assoc_with_invalid_allele_string():
    chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col = 'c', 'pos', 'pval', 'maf', 'rs', 'a1', 'a2'
    example = pd.DataFrame([
            [1, 'rs001', 100, 'G', 'T_INVALID', 0.01, 0.01],
            [1, 'rs002', 200, 'G', 'T', 0.01, 0.01]
        ]
        , columns = [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]
    )

    example_buff = io.StringIO()
    example.to_csv(example_buff, sep = '\t', header = True, index = False)
    example_buff.seek(0)
    logger = get_test_logger()
    iterable = plotpeaks.read_assoc(example_buff, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, logger)
    data = next(iterable)

    expected = pd.DataFrame([
            [1, 'rs002', 200, 'G', 'T', 0.01, 0.01]
        ]
        , columns = [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]
    )

    pd.testing.assert_frame_equal(expected, data)


def test_get_signals_when_no_signif_signals():
    """
    All detected signal's p-value is greater than significance threshold.
    Prevent error from occurring and just output an empty dataframe when this happens
    """
    chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col = 'c', 'pos', 'pval', 'maf', 'rs', 'a1', 'a2'
    example = pd.DataFrame([
            [1, 'rs001', 100, 'G', 'T', 0.01, '0.01'],
            [1, 'rs002', 200, 'G', 'T', 0.01, '0.01']
        ]
        , columns = [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]
    )

    example_buff = io.StringIO()
    example.to_csv(example_buff, sep = '\t', header = True, index = False)
    example_buff.seek(0)
    logger = get_test_logger()
    iterable = plotpeaks.read_assoc(example_buff, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, logger)


    signals = plotpeaks.get_signals(iterable, 5e-8, chr_col, pos_col, pval_col)
    assert signals.empty, 'Signals dataframe is not empty'


def test_read_asoc_for_NA_values_in_pval_column():
    """
    summary statistics file with NA values
    This is to check whether 'NA' value in the P-value column causes the read_assoc function to fail in reading the file.
    Return dataframe with row omitted.
    """
    chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col = 'c', 'pos', 'pval', 'maf', 'rs', 'a1', 'a2'
    example = pd.DataFrame([
            [1, 'rs001', 100, 'G', 'T', 0.01, 'NA'],
            [1, 'rs002', 200, 'G', 'T', 0.01, 0.01]
        ]
        , columns = [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]
    )
    example_buff = io.StringIO()
    example.to_csv(example_buff, sep = '\t', header = True, index = False)
    example_buff.seek(0)
    logger = get_test_logger()
    iterable = plotpeaks.read_assoc(example_buff, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, logger)
    data = next(iterable)
    expected = pd.DataFrame([
            [1, 'rs002', 200, 'G', 'T', 0.01, 0.01]
        ]
        , columns = [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]
    )

    pd.testing.assert_frame_equal(expected, data)


def test_read_asoc_for_extra_columns():
    """
    summary statistics file with unused column
    This is to check whether unused column is removed when read_assoc function reads in a summary stats file.
    Returns dataframe with omitted column.
    """
    chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, extra = 'c', 'pos', 'pval', 'maf', 'rs', 'a1', 'a2','extra'
    example = pd.DataFrame([
            [1, 'rs001', 100, 'G', 'T', 0.01, 0.01, 'NA'],
            [1, 'rs002', 200, 'G', 'T', 0.01, 0.01, 1]
        ]
        , columns = [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col, extra]
    )
    example_buff = io.StringIO()
    example.to_csv(example_buff, sep = '\t', header = True, index = False)
    example_buff.seek(0)
    logger = get_test_logger()
    iterable = plotpeaks.read_assoc(example_buff, chr_col, pos_col, pval_col, maf_col, rs_col, a1_col, a2_col, logger)
    data = next(iterable)
    expected = pd.DataFrame([
            [1, 'rs001', 100, 'G', 'T', 0.01, 0.01],
            [1, 'rs002', 200, 'G', 'T', 0.01, 0.01]
        ]
        , columns = [chr_col, rs_col, pos_col, a1_col, a2_col, maf_col, pval_col]
    )

    assert len(expected.columns) < len(data.columns)
