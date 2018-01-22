#!/usr/bin/env python
# -*- coding: utf-8 -*- 

#/nfs/team144/agilly/inter_manh_env/bin/python
import sys
from bokeh.plotting import *
from bokeh.models import *
import pandas as pd
import numpy as np
import logging
import json, requests, asr
from urllib.request import urlopen
from bokeh.io import output_file, show
from bokeh.layouts import gridplot
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
ff.columns=['ensembl_consequence', 'ensembl_rs', 'ps', 'ensembl_assoc']
ff['ensembl_assoc'].fillna("none", inplace=True)
ff=ff.groupby(ff['ps']).apply(lambda x: pd.Series({'ensembl_consequence' : ";".join(x['ensembl_consequence'].unique()), 'ensembl_rs' : ";".join(x['ensembl_rs'].unique()), 'ensembl_assoc' : ";".join(set(x['ensembl_assoc']))})).reset_index()


# Merge dataframes, drop signals that are not present in the dataset
emax=pd.merge(e, ff, on='ps', how='outer')
emax.loc[pd.isnull(emax['ensembl_rs']), 'ensembl_rs']="novel"
emax.loc[pd.isnull(emax['ensembl_consequence']), 'ensembl_consequence']="novel"
emax.dropna(subset=['chr'], inplace=True)
e=emax
e['ensembl_assoc'].fillna("none", inplace=True)

e['chr']=e['chr'].astype(int)


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





e=ColumnDataSource(e)
server = "http://ensembl.org" if build=="b38" else "http://grch37.ensembl.org";
url = server+"/Homo_sapiens/Variation/Explore?v=@ensembl_rs"
taptool = p.select(type=TapTool)
taptool.callback = OpenURL(url=url)
p.circle(sys.argv[3], 'logp', line_width=2, source=e, size=9, fill_color='col', line_color="black",  line_alpha='col_assoc')
p.xaxis[0].formatter.use_scientific = False

p2=figure(width=1500, height=300, x_range=p.x_range, tools=['tap'])
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
q=gridplot([[p], [p2]])
save(q)

