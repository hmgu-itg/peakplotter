import sys
import json
from typing import Tuple, Union, List, Dict

import requests
import pandas as pd

from peakplotter.data import CENTROMERE_B37, CENTROMERE_B38


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


def get_rsid_in_region(chrom, start, end, server, logger):

    logger.debug(f"get_variants_in_region({chrom}, {start}, {end}, {server}")
    snps = get_variants_in_region(chrom, start, end, server)
    logger.debug(f"get_phenos_in_region({chrom}, {start}, {end}, {server}")
    pheno = get_phenos_in_region(chrom, start, end, server)
    
    resp = make_resp(snps, pheno)
    return resp


def query_vep(chrom: pd.Series, pos: pd.Series, a1: pd.Series, a2: pd.Series, server: str, logger) -> List[Dict]:
    chrom = chrom.astype(str)
    pos = pos.astype(int).astype(str)
    a1 = a1.astype(str)
    a2 = a2.astype(str)
    queries: pd.Series = chrom+" "+pos+" . "+a1+" "+a2+" . . ."
    data = json.dumps({'variants': queries.to_list()})
    ext = "/vep/homo_sapiens/region"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    r = requests.post(server+ext, headers=headers, data=data)
    if not r.ok:
        logger.error(data)
        r.raise_for_status()
    return r.json()

# TODO: Merge with get_csq_novel_variants function
def _get_csq_novel_variants(e: pd.DataFrame, chrcol: str, pscol: str, a1col: str, a2col: str, server: str, logger) -> pd.DataFrame:
    """
    This function assumes that the input DataFrame object `e` has the following columns:
      - ps
      - ensembl_rs
      - ld
    """
    copied_e = e.copy()
    copied_e.loc[(
        copied_e['ensembl_rs']=="novel")
        & (copied_e[a1col]==copied_e[a2col])
        ,
        'ensembl_consequence']='double allele'
    novelsnps=copied_e.loc[(copied_e['ensembl_rs']=="novel") & (copied_e['ld']>0.1) & (copied_e['ensembl_consequence']!='double allele'),]
    if novelsnps.empty:
        return copied_e
    jData = query_vep(novelsnps[chrcol], novelsnps[pscol], novelsnps[a1col], novelsnps[a2col], server, logger)

    csq = pd.DataFrame(jData)

    for _, row in csq.iterrows():
        copied_e.loc[copied_e['ps']==row['start'], 'ensembl_consequence'] = row['most_severe_consequence']

    copied_e['ensembl_consequence'].replace('_', ' ')
    return copied_e

# TODO: Merge with _get_csq_novel_variants function
def get_csq_novel_variants(e, chrcol, pscol, a1col, a2col, server, logger):
    copied_e = e.copy()
    copied_e.loc[(copied_e['ensembl_rs']=="novel") & (copied_e[a1col]==copied_e[a2col]),'ensembl_consequence']='double allele'
    novelsnps=copied_e.loc[(copied_e['ensembl_rs']=="novel") & (copied_e['ld']>0.1) & (copied_e['ensembl_consequence']!='double allele'),]
    if novelsnps.empty:
        return copied_e
    pd.options.mode.chained_assignment = None # Temporarily suppress the SettingWithCopyWarning message
    novelsnps['query']=novelsnps[chrcol].astype(str)+" "+novelsnps[pscol].astype(int).astype(str)+" . "+novelsnps[a1col].astype(str)+" "+novelsnps[a2col].astype(str)+" . . ."
    pd.options.mode.chained_assignment = 'warn' # Reactivate the SettingWithCopyWarning message
    request='{ "variants" : ["'+'", "'.join(novelsnps['query'])+'" ] }'
    ext = "/vep/homo_sapiens/region"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    logger.info("\t\t\t🌐   Querying Ensembl VEP (POST) :"+server+ext)
    r = requests.post(server+ext, headers=headers, data=request)

    if not r.ok:
        logger.error("headers :"+request)
        r.raise_for_status()
        sys.exit(1)
    
    jData = json.loads(r.text)
    csq=pd.DataFrame(jData)


    for _, row in csq.iterrows():
        copied_e.loc[copied_e['ps']==row['start'],'ensembl_consequence']=row['most_severe_consequence']

    copied_e['ensembl_consequence'].replace('_', ' ')
    return copied_e


def get_overlap_genes(chrom, start, end, server, logger) -> pd.DataFrame:
    url = f'{server}/overlap/region/human/{chrom}:{start}-{end}?feature=gene'
    logger.info("\t\t\t🌐   Querying Ensembl overlap (Genes, GET) :"+url)
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
