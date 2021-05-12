#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

import io
import sys
import json
import logging

import asr
import requests
import numpy as np
import pandas as pd
from urllib.request import urlopen
from bokeh.io import output_file, show
from bokeh.models import *
from bokeh.layouts import gridplot
from bokeh.plotting import *

import helper_functions
from helper_functions import *

logging.basicConfig()
file=sys.argv[1]
outfile=file+".html"

chrcol=sys.argv[6]
pscol=sys.argv[3]
a1col=sys.argv[7]
a2col=sys.argv[8]
build=sys.argv[9]

pd.options.mode.chained_assignment = None  
d=pd.read_csv(file, sep=",",index_col=False)
output_file(outfile)

chrom=set(d[chrcol]) # This is for later to check whether the region window overlaps a centromere.
assert (len(chrom)==1),print("The chromosome spans across multiple chromosomes, which is impossible.")

# Check whether the region window overlaps a centromere. Print the information about the presence of centromeric region later in the script.
chrom=str(chrom.pop())
if build=="b38":
    url="https://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/38/Modeled_regions_for_GRCh38.tsv"
    s=requests.get(url).content
    cen_list=pd.read_csv(io.StringIO(s.decode('utf-8')), sep='\t', header=(0))
    if cen_list.empty==False:
        cen_list.drop(['HET7'], inplace=True)
        cen_list.drop(cen_list.columns[3], axis=1, inplace=True)
        cen_list.columns = ['chr','start','end']
        cen_start=int(cen_list[cen_list['chr'] == chrom]['start'])
        cen_end=int(cen_list[cen_list['chr'] == chrom]['end'])
    else:
        info("Failed to obtain GRCh38 centromere coordinates")
else:
   url="wget -O- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz | gunzip | grep -v -e chrX -e chrY | grep cen | mergeBed -i - | sed 's/chr//'"
   import subprocess
   sp=subprocess.check_output(url, shell=True)
   #print(sp.decode('utf-8'))
   cen_list=pd.read_table(io.StringIO(sp.decode('utf-8')), sep='\t', header=None)
   #print(chrom)
   if cen_list.empty==False:
      cen_list.apply(pd.to_numeric)
      #print(cen_list)
      cen_list.columns = ['chr','start','end']
      #chrom=pd.to_numeric(chrom)
      #print(cen_list.dtypes)
      cen_start=cen_list[cen_list['chr'] == np.int64(chrom)]['start']
      cen_end=cen_list[cen_list['chr'] == np.int64(chrom)]['end']
      cen_start=cen_start[0]
      cen_end=cen_end[0]
   else:
      info("Failed to obtain GRCh37 centromere coordinates")
region_start=min(d[pscol])
region_end=max(d[pscol])

cen_coordinate=set(range(cen_start,cen_end))
region=set(range(region_start, region_end))
cen_overlap=region.intersection(cen_coordinate)



# this was below before:     ("name", "   @"+sys.argv[4]),
hover= HoverTool(tooltips = [
        ("==============", "=============="),
    ("name", "   @rs"),
        ("RS-id", "@ensembl_rs"),
        ("ld","@ld"),
        ("M.A.F", "@"+sys.argv[5]),
        ("overlaps gene ", "@gene"),
        ("consequence ", "@ensembl_consequence"),
        ("known associations", "@ensembl_assoc")
])
p = figure(tools=[hover, 'box_zoom,wheel_zoom,reset,tap'], width=1500)
#e=d.query('ps > 33400000 and ps<33900000')
#print(d)
d.replace(r'\s+', np.nan, regex=True)
if (pscol != "ps"):
   d['ps']=d[pscol] 


#d.rename(columns={pscol:'ps'}, inplace=True)
e=d

# LD colouring
e['col']=pd.cut(e['ld'], 9, labels=range(1,10))
from bokeh.palettes import Spectral10
collol = pd.DataFrame({'pal':Spectral10})
e['col']= np.asarray(collol.ix[e['col']]).flatten()

# Log P column
e['logp']=-np.log10(e[sys.argv[2]])

## We get the list of rsids and phenotype associations in the region
server = "https://rest.ensembl.org" if build=="b38" else "http://grch37.rest.ensembl.org";
helper_functions.server=server;
ff=get_rsid_in_region(str(e[chrcol][0]), str(e[pscol].min()), str(e[pscol].max()))
print(ff.head())
ff.columns=['ensembl_consequence', 'ensembl_rs', 'ps', 'ensembl_assoc']
ff['ensembl_assoc'].fillna("none", inplace=True)
ff=ff.groupby(ff['ps']).apply(lambda x: pd.Series({'ensembl_consequence' : ";".join(x['ensembl_consequence'].unique()), 'ensembl_rs' : ";".join(x['ensembl_rs'].unique()), 'ensembl_assoc' : ";".join(set(x['ensembl_assoc']))})).reset_index()

print(e.head())

# Merge dataframes, drop signals that are not present in the dataset
emax=pd.merge(e, ff, on='ps', how='outer')
emax.loc[pd.isnull(emax['ensembl_rs']), 'ensembl_rs']="novel"
emax.loc[pd.isnull(emax['ensembl_consequence']), 'ensembl_consequence']="novel"
emax.dropna(subset=[chrcol], inplace=True)
e=emax
e['ensembl_assoc'].fillna("none", inplace=True)

e[chrcol]=e[chrcol].astype(int)


# Create the alpha vector for associations in LD
e['col_assoc']=0
e.loc[(e['ensembl_assoc']!="none") & (e['ld']>0.1), 'col_assoc']=1
e['col_assoc']=e['col_assoc'].astype(int)


# ENSEMBL consequences for variants in LD that do not have rs-ids
e=get_csq_novel_variants(e, chrcol, pscol, a1col, a2col)


# Below gets the genes > d
url = server+'/overlap/region/human/'+str(e[chrcol][0])+':'+str(int(e[pscol].min()))+'-'+str(int(e[pscol].max()))+'?feature=gene;content-type=application/json'
info("\t\t\t🌐   Querying Ensembl overlap (Genes, GET) :"+url)
response = urlopen(url).read().decode('utf-8')
jData = json.loads(response)
d=pd.DataFrame(jData)
e['gene']=""
for index, row in d.iterrows():
    e.loc[(e['ps']>row['start']) & (e['ps']<row['end']), 'gene']=e.loc[(e['ps']>row['start']) & (e['ps']<row['end']), 'gene']+";"+row['external_name']
e['gene']=e['gene'].str.replace(r'^\;', '')

ff=ff.loc[ff['ensembl_assoc']!="none",]
for index, row in ff.iterrows():
    span = Span(location=row['ps'],  line_color="firebrick", dimension='height')
    label = Label(x=row['ps'], y=e['logp'].max()+0.2, text=row['ensembl_assoc'], angle=90, angle_units="deg",  text_align="right", text_color="firebrick", text_font_size="11pt", text_font_style="italic")
    p.add_layout(label)
    p.add_layout(span)




 
e.to_csv(file+".csv", index=False)

e=ColumnDataSource(e)
server = "http://ensembl.org" if build=="b38" else "http://grch37.ensembl.org";
url = server+"/Homo_sapiens/Variation/Explore?v=@ensembl_rs"
taptool = p.select(type=TapTool)
taptool.callback = OpenURL(url=url)
p.circle(sys.argv[3], 'logp', line_width=2, source=e, size=9, fill_color='col', line_color="black",  line_alpha='col_assoc')
p.xaxis[0].formatter.use_scientific = False

p2=figure(width=1500, height=300, x_range=p.x_range, tools=['tap'])

if d.empty==False: # if d is not empty,
    ys=np.random.rand(len(d['end']))
    d['y']=ys

    d['color']="cornflowerblue"
    d['color'].ix[d['biotype']=="protein_coding"]="goldenrod"

    d['sens']="<"
    d['sens'].ix[d['strand']>0]=">"
    d['name']=d['sens']+d['external_name']
    p2.segment(x0=d['start'], x1=d['end'], y0=ys, y1=ys, line_width=4, color=d['color'])

    labels=LabelSet(x='start', y='y', text='name', source=ColumnDataSource(d))
    p2.add_layout(labels)
    p2.xaxis.visible = False
else:
    info("\t\t\t🌐   No genes overlap this genomic region.")

if (len(cen_overlap)>0): # Add indication of the centromeric region in the plot
    perc_overlap=int((len(cen_overlap)/len(region))*100)
    info("\t\t\t    {0}% of the genomic region overlaps a centromere".format(perc_overlap))

    xs=min(cen_overlap)
    xe=max(cen_overlap)+1
    cen_median=max(cen_overlap)-int((max(cen_overlap)-min(cen_overlap))/2)
    p2.segment(x0=xs, x1=xe, y0=0.5, y1=0.5, line_width=100, color="grey")
    d=dict(x=[cen_median],y=[0.9],text=["Centromeric_region"])
    labels=LabelSet(x="x", y="y", text="text",  source=ColumnDataSource(d))
    p2.add_layout(labels)




q=gridplot([[p], [p2]])
save(q)
