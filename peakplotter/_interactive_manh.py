from typing import Tuple, Union

import requests
import pandas as pd

from peakplotter import helper
from peakplotter.data import CENTROMERE_B37, CENTROMERE_B38

def _query(url, headers = None):
    if headers is None:
        headers = dict()
    headers['Content-Type'] = 'application/json'
    r = requests.get(url, headers = headers)
    if not r.ok:
        r.raise_for_status()

    decoded = r.json()
    return decoded

def get_variants_in_region(chrom, start, end, server) -> pd.DataFrame:
    # REST API Request
    url = f'{server}/overlap/region/human/{chrom}:{start}-{end}?feature=variation'
    decoded = _query(url)
    
    # Process json data to DataFrame
    snps = pd.DataFrame(decoded)
    snps['location'] = snps.seq_region_name.map(str)+":"+snps.start.map(str)+"-"+snps.end.map(str)
    reordered_snps = snps[['source', 'assembly_name', 'id',
                             'seq_region_name', 'start', 'end', 'location',
                             'strand', 'alleles', 'feature_type',
                             'consequence_type', 'clinical_significance']
                         ]
    
    return reordered_snps


def get_phenos_in_region(chrom, start, end, server) -> pd.DataFrame:
    url = f'{server}/phenotype/region/homo_sapiens/{chrom}:{start}-{end}?feature_type=Variation'
    json_data = _query(url)
    
    pheno = _process_pheno_r(json_data)
    
    return pheno
    

def _process_pheno_r(json_data) -> pd.DataFrame:
    """
    
    Example
    -------
    >>> empty_list = list()
    >>> _process_pheno_r(empty_list)
    """
    if len(json_data) == 0:
        return pd.DataFrame(columns = ['id', 'location', 'pheno'])
    
    data = list()
    for i in json_data:
        identifier = i['id']
        location = i['phenotype_associations'][0]['location']
        pheno = list()
        for assoc in i['phenotype_associations']:
            if assoc['source'] != 'COSMIC':
                pheno.append(assoc['description'])
        
        pheno = ';'.join(pheno)
        data.append([identifier, location, pheno])
    
    data = pd.DataFrame(data, columns = ['id', 'location', 'pheno'])
    assert not any(data.duplicated()), 'Duplicated data found'
    return data


def make_resp(snps: pd.DataFrame, pheno: pd.DataFrame) -> pd.DataFrame:
    resp = pd.merge(snps, pheno, on='location', how='outer')
    resp.drop(
        axis = 1 ,
        errors = 'ignore',
        labels = [
            "alleles", "assembly_name", "clinical_significance",
            "feature_type", "end", "seq_region_name",
            "strand", "source", "location", 'id_y'
        ],
        inplace = True
    )

    resp.rename(columns = {
        'start': 'ps',
        'id': 'rs',
        'id_x': 'rs',
        'consequence_type': 'consequence'},
    inplace = True)
    resp.dropna(inplace=True, subset=["rs"])
    resp['pheno'].fillna('none', inplace = True)
    return resp[['rs', 'ps', 'consequence', 'pheno']]


def get_rsid_in_region(chrom, start, end, server):
    print(f"[DEBUG] get_variants_in_region({chrom}, {start}, {end}, {server}")
    snps = get_variants_in_region(chrom, start, end, server)
    print(f"[DEBUG] get_phenos_in_region({chrom}, {start}, {end}, {server}")
    pheno = get_phenos_in_region(chrom, start, end, server)
    
    resp = make_resp(snps, pheno)
    return resp


def get_overlap_genes(chrom, start, end, server) -> pd.DataFrame:
    url = f'{server}/overlap/region/human/{chrom}:{start}-{end}?feature=gene'
    helper.info("\t\t\t🌐   Querying Ensembl overlap (Genes, GET) :"+url)
    decoded = _query(url)
    
    df = pd.DataFrame(decoded).fillna('')
    if 'external_name' not in df.columns:
        # NOTE: get_overlap_genes(1, 104210023, 105209701, 'https://rest.ensembl.org')
        # The above region resulted in a dataframe with no 'external_name'.
        # Can't write a test for this case due to the way the functions are written
        # so we modify the function without a test.
        df['external_name'] = ''
    return df

def _get_build_centromere_file(build: Union[int, str]) -> str:
    mapper = {
            'b38': CENTROMERE_B38,
            '38': CENTROMERE_B38,
            38: CENTROMERE_B38,
            'b37': CENTROMERE_B37,
            '37': CENTROMERE_B37,
            37: CENTROMERE_B37
        }
    return mapper[build]


def get_centromere_region(chrom: int, build: Union[str, int]) -> Tuple[int, int]:
    centromere = _get_build_centromere_file(build)

    start = centromere.loc[centromere['chrom'] == chrom, 'start'].to_list()[0]
    end = centromere.loc[centromere['chrom'] == chrom, 'end'].to_list()[0]

    return start, end
