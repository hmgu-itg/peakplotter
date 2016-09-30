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
import urllib2
from bokeh.io import gridplot, output_file, show

logging.basicConfig()
file=sys.argv[1]
outfile=file+".html"


pd.options.mode.chained_assignment = None  
d=pd.read_csv(file, sep=",",index_col=False)
output_file(outfile)
hover= HoverTool(tooltips = [
        ("==============", "=============="),
    ("name", "   @"+sys.argv[4]),
        ("RS-id", "@id"),
        ("ld","@ld"),
        ("M.A.F", "@"+sys.argv[5]),
        ("overlaps gene ", "@gene"),
        ("consequence ", "@consequence"),
        ("known associations", "@assoc")
])
p = figure(tools=[hover, 'box_zoom,wheel_zoom,reset,tap'], width=1500)
#e=d.query('ps > 33400000 and ps<33900000')
#print(d)
#d.replace(r'\s+', np.nan, regex=True)

e=d
#print(e)
e['col']=pd.cut(e['ld'], 9, labels=range(1,10))
from bokeh.palettes import Spectral10
collol = pd.DataFrame({'pal':Spectral10})

e['assoc']=e['assoc'].astype(basestring)


e['col_assoc']=e['assoc'].astype(basestring)

#e['col_assoc'].astype(basestring)


#e.replace(np.nan,0, regex=True)


e['col_assoc'].ix[pd.isnull(e['col_assoc'])]=0


e['col_assoc'].ix[e['col_assoc']=="none"]=0



e['col_assoc'].ix[e['col_assoc']!=0]=1
#print(e['col_assoc'])
e['col']= np.asarray(collol.ix[e['col']]).flatten()
e['logp']=-np.log10(e[sys.argv[2]])

print("Query1:Genes")
url = 'http://rest.ensembl.org/overlap/region/human/'+str(e['chr'][0])+':'+str(e['ps'].min())+'-'+str(e['ps'].max())+'?feature=gene;content-type=application/json'
print("Querying Ensembl with region "+url)
response = urllib2.urlopen(url).read()
jData = json.loads(response)
d=pd.DataFrame(jData)


print("Query2:GWASCAT")
url = 'http://rest.ensembl.org/overlap/region/human/'+str(e['chr'][0])+':'+str(e['ps'].min())+'-'+str(e['ps'].max())+'?feature=variation;content-type=application/json;variant_set=ph_nhgri'
print("Querying Ensembl with region "+url)
response = urllib2.urlopen(url).read()
jData = json.loads(response)
cat=pd.DataFrame(jData)

print("Query2:GWASCAT (POST)")

header=pd.DataFrame()
header['ids']=cat['id']

data = header.to_dict(orient="list")




server = "https://rest.ensembl.org"
ext = "/variation/homo_sapiens?phenotypes=1"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
r = requests.post(server+ext, headers=headers, data=json.dumps(data))
decoded = r.json()
s=repr(decoded)
#print(s)
for rsid in decoded:
	p.ray(x=decoded[rsid]['mappings'][0]['end'], color="firebrick", length=0, angle=90, angle_units="deg")
	label = Label(x=decoded[rsid]['mappings'][0]['end'], y=e['logp'].max(), text=decoded[rsid]['phenotypes'][0]['trait'], angle=90, angle_units="deg", text_align="right", text_color="firebrick", text_font_size="11", text_font_style="italic")
	p.add_layout(label)
#df=pd.DataFrame(decoded).transpose()



e=ColumnDataSource(e)

url = "http://ensembl.org/Homo_sapiens/Variation/Explore?v=@id"
taptool = p.select(type=TapTool)
taptool.callback = OpenURL(url=url)
p.circle(sys.argv[3], 'logp', line_width=2, source=e, size=9, fill_color='col', line_color="black",  line_alpha='col_assoc')


p2=figure(width=1500, height=300, x_range=p.x_range, tools=['tap'])
ys=np.random.rand(len(d['end']))
d['y']=ys
#d['center']=(d['start']+d['end'])/2

d['color']="cornflowerblue"
d['color'].ix[d['biotype']=="protein_coding"]="goldenrod"

d['sens']="<"
d['sens'].ix[d['strand']>0]=">"
d['name']=d['sens']+d['external_name']
p2.segment(x0=d['start'], x1=d['end'], y0=ys, y1=ys, line_width=4, color=d['color'])

# url = "http://ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=@id"
# taptool2 = p2.select(type=TapTool)
# taptool2.callback = OpenURL(url=url)
#p2.add_layout(Arrow(end=VeeHead(size=10), x_start=d['start'].astype(float), x_end=d['end'].astype(float), y_start=ys, y_end=ys, line_width=4))
#p2.add_layout(Arrow(end=VeeHead(size=10), x_start=15917698, x_end=15917698+10000, y_start=0, y_end=0, line_width=4))
labels=LabelSet(x='start', y='y', text='name', source=ColumnDataSource(d))
p2.add_layout(labels)
q=gridplot([[p], [p2]])
save(q)

