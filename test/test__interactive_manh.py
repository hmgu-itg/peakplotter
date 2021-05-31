
import numpy as np
import pandas as pd
from pandas import testing

from peakplotter._interactive_manh import make_resp, get_centromere_region

def test_make_resp_when_pheno_df_is_empty():

    example_snps = pd.DataFrame([
        ['dbSNP', 'GRCh38', 'rs100',
        '22', 200, 200,
        '22:200-200', 1, ['A', 'T'],
        'variation', 'intergenic_variant', []]
    ],
        columns = ['source', 'assembly_name', 'id', 'seq_region_name', 'start', 'end',
    'location', 'strand', 'alleles', 'feature_type', 'consequence_type',
    'clinical_significance']
    )

    example_empty_phenos = pd.DataFrame(
        columns = ['phenotype_associations', 'id', 'pheno', 'location']
        )

    expected = pd.DataFrame([
        ['rs100', 200, 'intergenic_variant', 'none']],
        columns = ['rs', 'ps', 'consequence', 'pheno']
    )


    testing.assert_frame_equal(expected, make_resp(example_snps, example_empty_phenos))



def test_make_resp_when_phenotype_present():
    example_snps = pd.DataFrame([
        ['dbSNP', 'GRCh38', 'rs100',
        '22', 100, 100,
        '22:100-100', 1, list(['G', 'A']),
        'variation', 'synonymous_variant', list(['likely benign'])]
        ], 
        columns = ['source', 'assembly_name', 'id', 'seq_region_name', 'start', 'end',
        'location', 'strand', 'alleles', 'feature_type', 'consequence_type',
        'clinical_significance'])

    example_processed_phenos = pd.DataFrame([
                [list([{'attributes': {'associated_gene': 'DGCR2'}, 'source': 'HGMD-PUBLIC', 'description': 'Annotated by HGMD', 'location': '22:19041168-19041168'}]),
                'CM119374',
                'Annotated by HGMD',
                '22:100-100']],
                columns = ['phenotype_associations', 'id', 'pheno', 'location'])

    expected = pd.DataFrame([
        ['rs100', 100, 'synonymous_variant', 'Annotated by HGMD'],
        ], columns = ['rs', 'ps', 'consequence', 'pheno']
        )

    testing.assert_frame_equal(expected, make_resp(example_snps, example_processed_phenos))


def test_get_centromere_region_build():
    start, end = get_centromere_region(1, build = 38)
    assert start == 122026460
    assert end == 125184587
    
    start, end = get_centromere_region(22, build = 38)
    assert start == 12954789
    assert end == 15054318


    start, end = get_centromere_region(1, build = 37)
    assert start == 121500000
    assert end == 128900000
    
    start, end = get_centromere_region(22, build = 37)
    assert start == 12200000
    assert end == 17900000
