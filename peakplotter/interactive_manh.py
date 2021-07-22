#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

import sys

import numpy as np
import pandas as pd
from bokeh.io import output_file
from bokeh.models import HoverTool, TapTool, OpenURL, LabelSet, Span, Label, ColumnDataSource
from bokeh.layouts import gridplot
from bokeh.plotting import figure, save
from bokeh.palettes import Spectral10

from peakplotter import helper # TODO: Change this back to relative import after we replace plotpeaks.sh to a python equivalent.
from peakplotter import _interactive_manh

def interactive_manh(file, pvalcol, pscol, mafcol, chrcol, a1col, a2col, build: str, logger):
    outfile=file+".html"

    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)
    pd.options.mode.chained_assignment = None  
    d=pd.read_csv(file, sep=",",index_col=False)
    output_file(outfile)

    chrom=set(d[chrcol]) # This is for later to check whether the region window overlaps a centromere.
    assert len(chrom)==1, "The chromosome spans across multiple chromosomes, which is impossible."

    # Check whether the region window overlaps a centromere. 
    # Print the information about the presence of centromeric region later in the script.
    chrom = int(chrom.pop())

    cen_start, cen_end = _interactive_manh.get_centromere_region(chrom, build)
    
    region_start = min(d[pscol])
    region_end = max(d[pscol])

    cen_coordinate = set(range(cen_start, cen_end))
    region = set(range(region_start, region_end))
    cen_overlap = region.intersection(cen_coordinate)



    # this was below before:     ("name", "   @"+sys.argv[4]),
    hover = HoverTool(tooltips = [
            ("==============", "=============="),
            ("name", "   @rs"),
            ("RS-id", "@ensembl_rs"),
            ("ld","@ld"),
            ("M.A.F", "@"+mafcol),
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

    collol = pd.DataFrame({'pal':Spectral10})
    # e['col']= np.asarray(collol.ix[e['col']]).flatten() # .ix deprecated
    e['col']= np.asarray(collol.loc[e['col']]).flatten()

    # Log P column
    e['logp']=-np.log10(e[pvalcol])

    start = int(e[pscol].min())
    end = int(e[pscol].max())
    ## We get the list of rsids and phenotype associations in the region
    server = _interactive_manh.get_build_server(build)
    logger.debug(f"get_rsid_in_region({chrom}, {start}, {end}, {server}, logger")
    ff = _interactive_manh.get_rsid_in_region(chrom, start, end, server, logger)
    logger.debug(ff.head())
    columns = ['ensembl_rs', 'ps', 'ensembl_consequence', 'ensembl_assoc']
    ff.columns = columns

    dedup_ff = ff[~ff['ps'].duplicated(keep = False)]
    dup_ff = ff[ff['ps'].duplicated(keep = False)]
    grouped_dup_ff = dup_ff.groupby('ps').apply(lambda x: pd.Series({
                        'ensembl_consequence': ";".join(x['ensembl_consequence'].unique()),
                        'ensembl_rs': ";".join(x['ensembl_rs'].unique()),
                        'ensembl_assoc': ";".join(set(x['ensembl_assoc']))}))\
                    .reset_index()

    grouped_ff = pd.concat([
                    dedup_ff[columns],
                    grouped_dup_ff[columns]
                ]).sort_values('ps').reset_index(drop = True)
    logger.debug(e.head())

    # Merge dataframes, drop signals that are not present in the dataset
    emax=pd.merge(e, grouped_ff, on='ps', how='outer')
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
    logger.debug(f"e=helper.get_csq_novel_variants(e, '{chrcol}', '{pscol}', '{a1col}', '{a2col}', '{server}')")
    e = helper.get_csq_novel_variants(e, chrcol, pscol, a1col, a2col, server)


    # Below gets the genes > d
    
    d = _interactive_manh.get_overlap_genes(chrom, start, end, server, logger)
    e['gene']=""
    for index, row in d.iterrows():
        external_name = row['external_name']
        if external_name=='':
            # TODO: This case is things like pseudogenes, IncRNA, and rRNA genes.
            # Check with Arthur whether these genes should just be removed completely
            external_name = row['biotype']
        e.loc[(e['ps']>row['start']) & (e['ps']<row['end']), 'gene']=e.loc[(e['ps']>row['start']) & (e['ps']<row['end']), 'gene']+";"+external_name
    e['gene']=e['gene'].str.replace(r'^\;', '')

    grouped_ff=grouped_ff.loc[grouped_ff['ensembl_assoc']!="none",] # TODO: Currently, the .fillna('none') is done inside the _interactive_manh.make_resp function. This would be very difficult to understand for others. Change this to something simpler.
    # for index, row in grouped_ff.iterrows():
        # Make red line with previously associated information
        # TODO: Re-enable this feature later possibly with different modes (all, only functional, none),
        # or by allowing the viewer interactively enable and disable their choice of traits within the bokeh plot.
        # https://docs.bokeh.org/en/latest/docs/user_guide/interaction/legends.html
        # https://docs.bokeh.org/en/latest/docs/user_guide/interaction/widgets.html
        # span = Span(location=row['ps'],  line_color="firebrick", dimension='height')
        # label = Label(x=row['ps'], y=e['logp'].max()+0.2, text=row['ensembl_assoc'], angle=90, angle_units="deg",  text_align="right", text_color="firebrick", text_font_size="11pt", text_font_style="italic")
        # p.add_layout(label)
        # p.add_layout(span)

    
    e.to_csv(file+".csv", index=False)

    e=ColumnDataSource(e)
    server = "http://ensembl.org" if build=="b38" else "http://grch37.ensembl.org"
    url = f"{server}/Homo_sapiens/Variation/Explore?v=@ensembl_rs"
    taptool = p.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    p.circle(pscol, 'logp', line_width=2, source=e, size=9, fill_color='col', line_color="black",  line_alpha='col_assoc')
    p.xaxis[0].formatter.use_scientific = False

    # Make bottom figure with genes
    p2=figure(width=1500, height=300, x_range=p.x_range, tools=['tap'])
    if d.empty==False: # if d is not empty,
        ys=np.random.rand(len(d['end']))
        d['y']=ys

        d['color']="cornflowerblue"
        # d['color'].ix[d['biotype']=="protein_coding"]="goldenrod"
        d['color'].loc[d['biotype']=="protein_coding"]="goldenrod"

        d['sens']="<"
        # d['sens'].ix[d['strand']>0]=">"
        d['sens'].loc[d['strand']>0]=">"
        d['name']=d['sens']+d['external_name']
        p2.segment(x0=d['start'], x1=d['end'], y0=ys, y1=ys, line_width=4, color=d['color'])

        labels=LabelSet(x='start', y='y', text='name', source=ColumnDataSource(d))
        p2.add_layout(labels)
        p2.xaxis.visible = False
    else:
        logger.info("\t\t\t🌐   No genes overlap this genomic region.")

    # TODO: Move all the centromere code down here.
    if (len(cen_overlap)>0): # Add indication of the centromeric region in the plot
        perc_overlap=int((len(cen_overlap)/len(region))*100)
        logger.info("\t\t\t    {0}% of the genomic region overlaps a centromere".format(perc_overlap))

        xs=min(cen_overlap)
        xe=max(cen_overlap)+1
        cen_median=max(cen_overlap)-int((max(cen_overlap)-min(cen_overlap))/2)
        p2.segment(x0=xs, x1=xe, y0=0.5, y1=0.5, line_width=100, color="grey")
        d=dict(x=[cen_median],y=[0.9],text=["Centromeric_region"])
        labels=LabelSet(x="x", y="y", text="text",  source=ColumnDataSource(d))
        p2.add_layout(labels)

    q=gridplot([[p], [p2]])
    save(q)


if __name__ == '__main__':
    file=sys.argv[1]
    pvalcol=sys.argv[2]
    pscol=sys.argv[3]
    rscol=sys.argv[4]
    mafcol=sys.argv[5]
    chrcol=sys.argv[6]
    a1col=sys.argv[7]
    a2col=sys.argv[8]
    build=sys.argv[9]
    
    interactive_manh(file, pvalcol, pscol, rscol, mafcol, chrcol, a1col, a2col, build)
