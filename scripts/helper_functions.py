#!/usr/bin/env python3

import re
import sys
import json
import subprocess
import urllib.request
import urllib.parse

import asr
import requests
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# number by which to multiply the normal number of requests to Ensembl
# Increase if lots of 504 Timeout errors
ENSEMBL_USELESSNESS_COEFFICIENT=3

def info (*strs):
	outstr="[INFO]";
	for string in strs:
		outstr+=" "+str(string)
	print(outstr);
	return;

def warn (*strs):
	outstr="[WARNING]";
	for string in strs:
		outstr+=" "+str(string)
	print(outstr, file=sys.stderr);
	return;

class GeneCoordinates:
	chrom=0
	start=0
	end=0
	gene_id=""
	name=""
	def __init__(self, chrom, start, end, gene_id, name):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.gene_id=gene_id
		self.name=name

	def extend(self, margin):
		self.start-=int(margin);
		self.end+=int(margin);


def get_csq_novel_variants(e, chrcol, pscol, a1col, a2col):
	global server
	e.loc[(e['ensembl_rs']=="novel") & (e[a1col]==e[a2col]),'ensembl_consequence']='double allele'
	novelsnps=e.loc[(e['ensembl_rs']=="novel") & (e['ld']>0.1) & (e['ensembl_consequence']!='double allele'),]
	if novelsnps.empty:
		return e
	novelsnps['query']=novelsnps[chrcol].astype(str)+" "+novelsnps[pscol].astype(int).astype(str)+" . "+novelsnps[a1col].astype(str)+" "+novelsnps[a2col].astype(str)+" . . ."
	request='{ "variants" : ["'+'", "'.join(novelsnps['query'])+'" ] }'
	ext = "/vep/homo_sapiens/region"
	headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	info("\t\t\t🌐   Querying Ensembl VEP (POST) :"+server+ext)
	r = requests.post(server+ext, headers=headers, data=request)
	 
	if not r.ok:
		print("headers :"+request)
		r.raise_for_status()
		sys.exit()
	 
	jData = json.loads(r.text)
	csq=pd.DataFrame(jData)


	for index,row in csq.iterrows():
	    e.loc[e['ps']==row['start'],'ensembl_consequence']=row['most_severe_consequence']

	e['ensembl_consequence'].replace('_', ' ')
	return e

def get_coordinates(gene_name):
    '''
    Function to return the genomic coordinates of a gene submitted (GRCh37)
    output data: chromosome, start, end, stable ID
    '''

    global server
    ext = "/lookup/symbol/homo_sapiens/%s?expand=0" % (gene_name)

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    gc=GeneCoordinates(int(decoded["seq_region_name"]), int(decoded["start"]), int(decoded["end"]), decoded["id"], gene_name)
    return(gc)

def get_rsid_in_region(c, start, end):
	url = server+'/overlap/region/human/'+str(c)+":"+str(start)+"-"+str(end)+'?feature=variation;content-type=application/json;';
	info("\t\t\t🌐   Querying Ensembl with region "+url)
	response = urllib.request.urlopen(url).read().decode('utf-8')
	jData = json.loads(response)
	snps=pd.DataFrame(jData)
	snps['location']=snps.seq_region_name.map(str)+":"+snps.start.map(str)+"-"+snps.end.map(str)

	url = server+'/phenotype/region/homo_sapiens/'+str(c)+":"+str(start)+"-"+str(end)+'?feature_type=Variation;content-type=application/json;';
	info("\t\t\t🌐   Querying Ensembl with region "+url)
	response = urllib.request.urlopen(url).read().decode('utf-8')
	jData = json.loads(response)
	pheno=pd.DataFrame(jData)
	#print(pheno['phenotype_associations'])
	pheno['pheno']=""
	pheno['location']=""
	for index, variant in pheno.iterrows():
		for assoc in variant.phenotype_associations:
			if assoc['source'] != 'COSMIC':
				try:
					variant.pheno=";".join(set(variant.pheno.split(";").append(assoc['description'])))
				except TypeError as e:
					variant.pheno=assoc['description']
				
				#variant.pheno=assoc['description'] if (variant.pheno=="") else variant.pheno+";"+assoc['description'];
				
				variant.location=assoc['location']
	resp=pd.merge(snps, pheno, on='location', how='outer')
	if pheno.empty==True:
		resp.drop(["alleles", "assembly_name", "clinical_significance", "feature_type", "end", "seq_region_name", "strand", "source", "location"], axis=1, inplace=True)
		resp.rename(columns = {'start':'ps', 'id':'rs', 'consequence_type':'consequence'}, inplace = True)	
	else:
		resp.drop(["alleles", "assembly_name", "clinical_significance", "feature_type", "end", "seq_region_name", "phenotype_associations", "strand", "source", "id_y", "location"], axis=1, inplace=True)
		resp.rename(columns = {'start':'ps', 'id_x':'rs', 'consequence_type':'consequence'}, inplace = True)	
	resp.dropna(inplace=True, subset=["rs"])
	return(resp)

def get_rsid_in_region_old(gc):
	c=gc.chrom
	start=gc.start
	end=gc.end
	url = server+'/overlap/region/human/'+str(c)+":"+str(start)+"-"+str(end)+'?feature=variation;content-type=application/json;';
	info("Querying Ensembl with region "+url)
	response = urllib.request.urlopen(url).read().decode('utf-8')
	jData = json.loads(response)
	cat=pd.DataFrame(jData)
	nsplit=ENSEMBL_USELESSNESS_COEFFICIENT*int(cat.count()[0]/1000)
	resp=pd.DataFrame(columns=('rs', 'ps', 'consequence', 'pheno'))
	info("\t\t\t🌐   Performing "+str(nsplit)+" phenotype requests...");
	ext = "/variation/homo_sapiens?phenotypes=1"
	j=0
	for i in np.array_split(cat['id'], nsplit):
		j=j+1
		info("\t\t\t⌛   GWASCAT (POST) at "+server+ext+" ("+str(j)+'/'+str(nsplit)+")")
		header=pd.DataFrame()
		header['ids']=i
		data = header.to_dict(orient="list")
		headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
		data=json.dumps(data)
		data = data.encode('ascii')
		req = urllib.request.Request(server+ext, data=data, headers=headers)
		with urllib.request.urlopen(req) as response:
			r=response.read().decode('utf-8')
			jData = json.loads(r)
			#pp.pprint(jData)
			for rsid in jData:
				ps=jData[rsid]['mappings'][0]['end']
				consequence=jData[rsid]['most_severe_consequence']
				pheno=""
				for k in jData[rsid]['phenotypes']:
					pheno=k['trait'] if (pheno=="") else pheno+";"+k['trait'];
				resp=resp.append({'rs':rsid, 'ps':int(ps), 'consequence':consequence, 'pheno':pheno}, ignore_index=True)
	return(resp)

def fetch_single_point(gc, sp_results):
	c=gc.chrom
	start=gc.start
	end=gc.end
	sp = pd.DataFrame();
	task = subprocess.Popen(["tabix", "-h", sp_results, str(c)+":"+str(start)+"-"+str(end)], stdout=subprocess.PIPE);
	sp=pd.read_table(task.stdout);
	task = subprocess.Popen(["zgrep", "-m1", "chr", sp_results], stdout=subprocess.PIPE);
	sp.columns=task.stdout.read().decode('UTF-8').split();
	return(sp)

def read_variants_from_gene_set(gc, input_monster):
	c=gc.chrom
	variants=pd.read_table(input_monster)
	variants.columns=[w.replace('chr'+str(c)+"_", '') for w in variants.columns]
	variants.columns=[re.sub("_.*", '', w) for w in variants.columns]
	variants.drop(variants.columns[[0,1]], axis=1, inplace=True)
	variants=variants.transpose()
	variants['ps']=variants.index
	variants.index=range(variants.count()[0])
	if (len(variants.columns)==2):
		variants.columns=["weight", "ps"]
	else:
		#Sometimes runs have no weights
		variants.columns=["ps"]
		variants['weight']=1
	variants.ps=pd.to_numeric(variants.ps, errors='coerce')
	return(variants)
