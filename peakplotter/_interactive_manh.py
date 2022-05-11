import json
import itertools
from typing import Tuple, Union

import requests
import pandas as pd

from peakplotter.data import CENTROMERE_B37, CENTROMERE_B38


_ENSEMBL_MAX = 5_000_000
_ENSEMBL_PARTS = 2_000_000
# Ensembl REST API GET for some region based queries limit the
# query range to 5Mb. 

def get_build_server(build: Union[int, str]) -> str:
	B38_SERVER = "https://rest.ensembl.org"
	B37_SERVER = "http://grch37.rest.ensembl.org"
	mapper = {
			'b38': B38_SERVER,
			'38': B38_SERVER,
			38: B38_SERVER,
			'b37': B37_SERVER,
			'37': B37_SERVER,
			37: B37_SERVER
		}
	return mapper[build]


def _query(url, headers = None):
    if headers is None:
        headers = dict()
    headers['Content-Type'] = 'application/json'
    r = requests.get(url, headers = headers)
    if not r.ok:
        r.raise_for_status()

    decoded = r.json()
    return decoded


def _divide_query_parts(start: int, end: int) -> list:
    remain = end - start
    pos = start
    parts = list()
    while remain:
        if remain // _ENSEMBL_MAX:
            parts.append((pos, pos+_ENSEMBL_PARTS))
            pos+=_ENSEMBL_PARTS
            remain-=_ENSEMBL_PARTS
        else:
            parts.append((pos, pos+remain))
            break
    return parts


def get_variants_in_region(chrom, start, end, server) -> pd.DataFrame:
    # REST API Request
    if end - start > _ENSEMBL_MAX:
        parts = _divide_query_parts(start, end)
        decoded = list()
        for (start, end) in parts:
            url = f'{server}/overlap/region/human/{chrom}:{start}-{end}?feature=variation'
            decoded.extend(_query(url))
    else:
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
    if end - start > _ENSEMBL_MAX:
        parts = _divide_query_parts(start, end)
        json_data = list()
        for (start, end) in parts:
            url = f'{server}/phenotype/region/homo_sapiens/{chrom}:{start}-{end}?feature_type=Variation'
            json_data.extend(_query(url))
    else:
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


def get_rsid_in_region(chrom, start, end, server, logger):

    logger.debug(f"get_variants_in_region({chrom}, {start}, {end}, '{server}')")
    snps = get_variants_in_region(chrom, start, end, server)
    logger.debug(f"get_phenos_in_region({chrom}, {start}, {end}, '{server}')")
    pheno = get_phenos_in_region(chrom, start, end, server)
    
    resp = make_resp(snps, pheno)
    return resp


def get_csq(data: pd.DataFrame, build: int = 38) -> pd.DataFrame:
    '''
    Queries Ensembl VEP for the input variants.
    
    Parameters
    ----------
    data : pd.DataFrame
        Expects the dataframe to have the following columns: chrom, pos, a1, and a2.
    build : int
        GRCh build number to query. Choices are 37 and 38

    Returns
    -------
    pd.DataFrame
    '''

    server = get_build_server(build)
    url = server + "/vep/homo_sapiens/region"

    headers = {"Content-Type" : "application/json", "Accept": "application/json"}

    # Query VEP in chunks of 200 variants (because Ensembl limits POST sizes to 200)
    vep_data = dict()
    for i in range(0, data.shape[0], 200):
        subset = data.iloc[i:i+200]
        json_data = {'variants': [f'{row[0]} {row[1]} . {row[2]} {row[3]} . . .' for _, row in subset.iterrows()]}
        json_data = json.dumps(json_data)
        r = requests.post(url, headers=headers, data=json_data)
        if not r.ok:
            r.raise_for_status()
        decoded = {i['input']: i for i in r.json()}
        vep_data.update(decoded)

    # Extract the consequence info from the query result
    _output = list()
    for k, v in vep_data.items():
        if 'transcript_consequences' in v:
            csqs = ';'.join(set(itertools.chain(*[i['consequence_terms'] for i in v['transcript_consequences']])))
        else:
            csqs = v['most_severe_consequence']
        _output.append([k, csqs])
    output = pd.DataFrame(_output, columns = ['input', 'csq'])
    inputs = output['input'].str.split(' ')
    output.insert(0, 'chrom', inputs.str[0])
    output.insert(1, 'ps', inputs.str[1].astype(int))
    output.insert(2, 'a1', inputs.str[3])
    output.insert(3, 'a2', inputs.str[4])
    output.drop(columns = 'input', inplace = True)
    return output




def get_overlap_genes(chrom, start, end, server) -> pd.DataFrame:
    if end - start > _ENSEMBL_MAX:
        parts = _divide_query_parts(start, end)
        decoded = list()
        for (start, end) in parts:
            url = f'{server}/overlap/region/human/{chrom}:{start}-{end}?feature=gene'
            decoded.extend(_query(url))
    else:
        url = f'{server}/overlap/region/human/{chrom}:{start}-{end}?feature=gene'
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
