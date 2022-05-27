#!/usr/bin/env python3

import functools
from datetime import datetime

import numpy as np
import pandas as pd
from bokeh.io import output_file
from bokeh.core.properties import Int
from bokeh.models import (HoverTool, TapTool, OpenURL, 
                         LabelSet, ColumnDataSource, Title, 
                         CustomJS, RangeSlider, CustomJSFilter, 
                         CDSView, CheckboxButtonGroup)
from bokeh.models.formatters import NumeralTickFormatter
from bokeh.layouts import gridplot
from bokeh.plotting import Figure, save
from bokeh.palettes import Spectral10

from . import __version__
from . import _interactive_manh, _ensembl_consequence


class GenomeView(Figure):
    __subtype__ = "GenomeView"
    __view_module__ = "bokeh" # https://github.com/bokeh/bokeh/issues/9412
    __view_model__ = "Plot"
    
    def __init__(self, *args, width = 1500, **kw):
        hover = HoverTool(tooltips = [
            ("==============", "=============="),
            ("RS-id", "@ensembl_rs"),
            ("pos", "@ps"),
            ("(a1, a2)", "@a1, @a2"),
            ("M.A.F", "@maf"),
            ("p-value", "@{p-value}"),
            ("ld", "@ld"),
            ("overlaps gene", "@gene"),
            ("consequence", "@ensembl_consequence"),
            ("known associations", "@ensembl_assoc")
        ])
        
        if 'tools' not in kw:
            kw['tools'] = [hover, 'box_zoom,wheel_zoom,undo,redo,reset,tap']
        
        super().__init__(*args, width = width, **kw)
        
        self.xaxis.formatter = NumeralTickFormatter()
        self.toolbar_location="above"
        
    def update_view(self, data, build = 38):
        source = ColumnDataSource(data)
        server = "http://ensembl.org" if build==38 else "http://grch37.ensembl.org"
        url = f"{server}/Homo_sapiens/Variation/Explore?v=@ensembl_rs"
        taptool = self.select(type=TapTool)
        taptool.callback = OpenURL(url=url)
        self.circle(
            'ps', 
            'logp', 
            line_width = 2, 
            source = source, 
            size = 9, 
            fill_color = 'col', 
            line_color = "black", 
            line_alpha = 'col_assoc')
        
        
        chrom = set(data['chrom']).pop()
        start = int(data['ps'].min())
        end = int(data['ps'].max())
        self.title = f'chr{chrom}:{start}-{end}'

        
class GeneView(Figure):
    __subtype__ = 'GeneView'
    __view_module__ = "bokeh"
    __view_model__ = "Plot"
    
    line_width = Int(4, help="The thickness of the gene plotted.")
    
    def __init__(self, *args, height = 500, width = 1500, **kw):
        if 'tools' not in kw:
            kw['tools'] = ['tap']
            
        super().__init__(*args, height = height, width = width, **kw)
        
        self.xaxis.formatter = NumeralTickFormatter()
        self.yaxis.visible = False

    def update_view(self, data):
        if self.renderers:
            for seg in self.renderers:
                seg.visible = False
            labels = self.select(LabelSet)
            for lab in labels:
                lab.visible = False
                
            
        genes = data.copy()
        
        gene_count = data.shape[0]
        rng = np.random.default_rng()
        ys = rng.choice(gene_count, size=gene_count, replace=False)

        genes['y'] = ys

        genes['color'] = "cornflowerblue"
        genes['protein_coding'] = 'not protein coding'
        genes.loc[genes['biotype']=="protein_coding", 'color']="goldenrod"
        genes.loc[genes['biotype']=="protein_coding", 'protein_coding']="protein coding"

        genes['sens'] = "<"
        genes.loc[genes['strand']>0, 'sens'] = ">"
        
        no_names = genes['external_name'] == ''
        genes.loc[no_names, 'external_name'] = genes.loc[no_names, 'gene_id']

        genes['name'] = genes['sens'] + genes['external_name']
        
        # source = ColumnDataSource(genes)
        
        for _type in ('not protein coding', 'protein coding'):
            subset_source = ColumnDataSource(genes.loc[genes['protein_coding']==_type])
            segments = self.segment(
                source = subset_source,
                x0='start',
                x1='end',
                y0='y',
                y1='y',
                line_width=self.line_width,
                color='color',
                legend_label = _type
            )
            
            labels = LabelSet(x='start', y='y', text='name', source = subset_source)
            self.add_layout(labels)
            
            segments.js_on_change('visible', CustomJS(args={'ls': labels}, code="ls.visible = cb_obj.visible;"))
            

        # labels = LabelSet(x='start', y='y', text='name', source=source)
        # self.add_layout(labels)
        self.legend.location = "top_left"
        self.legend.click_policy="hide"


def _make_grouped_ff(chrom, start, end, build, logger):
    server = _interactive_manh.get_build_server(build)
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
    
    grouped_ff['ps'] = grouped_ff['ps'].astype(int)
    
    return grouped_ff


def make_view_data(file, chrcol, pscol, a1col, a2col, pvalcol, mafcol, build, logger, vep_ld=0.1):
    d = pd.read_csv(file, sep=",", index_col=False)
    d.replace(r'\s+', np.nan, regex=True, inplace = True)
    d.rename(inplace = True, 
             columns = {
                 pvalcol: 'p-value',
                 chrcol: 'chrom',
                 pscol: 'ps',
                 a1col: 'a1',
                 a2col: 'a2',
                 mafcol: 'maf'
             }
    )

    chrom = set(d['chrom']) # This is for later to check whether the region window overlaps a centromere.
    assert len(chrom)==1, "The chromosome spans across multiple chromosomes, which is impossible."
    chrom = int(chrom.pop())
    start = int(d['ps'].min())
    end = int(d['ps'].max())
    server = _interactive_manh.get_build_server(build)

    # LD colour
    d['col'] = pd.cut(d['ld'], 9, labels = Spectral10[1:])

    d['logp'] = -np.log10(d['p-value'])
    
    grouped_ff = _make_grouped_ff(chrom, start, end, build, logger)


    d = pd.merge(d, grouped_ff, on='ps', how='left')
    d['ensembl_rs'].fillna('novel', inplace = True)
    d['ensembl_consequence'].fillna('novel', inplace = True)

    d['ensembl_assoc'].replace('', np.nan, inplace = True)
    d['ensembl_assoc'].fillna("none", inplace=True)
    # Remove duplicate entries of ensembl_assoc 
    d['ensembl_assoc'] = [';'.join(set(i.split(';')).difference({''})) for i in d['ensembl_assoc']]  
    d['ensembl_assoc'] = d['ensembl_assoc'].str.replace(r'^none;', '', regex=True)  
    # Create the alpha vector for associations in LD
    d['col_assoc']=0
    d.loc[d['ensembl_assoc']!="none", 'col_assoc']=1
    d['col_assoc'] = d['col_assoc'].astype(int)



    # ENSEMBL consequences for variants in LD that do not have rs-ids
    logger.info("\t\t\t🌐   Querying Ensembl VEP (POST)")
    d.loc[(d['ensembl_rs']=="novel") & (d['a1']==d['a2']), 'ensembl_consequence'] = 'double allele'

    if vep_ld<=0.0:
        vep_ld = -1.0
        # We do this because we sometimes assign LD=-0.01 to variants where plink can't calculate LD.
        # See line 273 of plotpeaks.py
    unknown_csqs = d.loc[(d['ensembl_consequence']=="novel") & (d['ld']>vep_ld) & (d['ensembl_consequence']!='double allele'), ['chrom', 'ps', 'a1', 'a2']]
    csqs = _interactive_manh.get_csq(unknown_csqs, build)
    csqs['chrom'] = csqs['chrom'].astype(d['chrom'].dtype)

    logger.debug(f'd.shape = {d.shape}')
    logger.debug(f'unknown_csqs.shape = {unknown_csqs.shape}')
    logger.debug(f'csqs.shape = {csqs.shape}')

    e = d.merge(csqs, how = 'left')
    e.loc[d['ensembl_consequence']=='novel', 'ensembl_consequence'] = e.loc[e['ensembl_consequence']=='novel', 'csq']
    e.drop(columns = 'csq', inplace = True)
    e['ensembl_consequence'].fillna('unknown(needs debugging)', inplace=True)
    
    e['ensembl_consequence_level'] = [min([_ensembl_consequence._consequences.get(i, 4) for i in v.split(';')]) for v in e['ensembl_consequence']]

    logger.debug(f'e.shape = {e.shape}')

    genes = _interactive_manh.get_overlap_genes(chrom, start, end, server)
    f_genes = genes[genes['external_name']!='']
    e['gene']=""
    for index, row in f_genes.iterrows():
        external_name = row['external_name']
        if external_name=='':
            # TODO: This case is things like pseudogenes, IncRNA, and rRNA genes.
            # Check with Arthur whether these genes should just be removed completely
            external_name = row['biotype']
        overlap = (e['ps']>row['start']) & (e['ps']<row['end'])
        e.loc[overlap, 'gene']=e.loc[overlap, 'gene']+";"+external_name
    e['gene'] = e['gene'].str.replace(r'^\;', '', regex=True)
    
    return e, genes


def _create_peakplot(e, genes, build, logger):
    ## Make GenomeView plot
    server = "http://grch37.ensembl.org" if build==37 or build=="b37" or build=="37" else "http://ensembl.org"
    ensembl_rs_url = f"{server}/Homo_sapiens/Variation/Explore?v=@ensembl_rs"
    genome = GenomeView()
    genome.select(type=TapTool)
    genome.callback = OpenURL(url=ensembl_rs_url)

    range_slider = RangeSlider(start = -0.01, end = 1.0, value = (-0.01, 1.0), step = 0.01, title = 'LD')

    checkbox_button = CheckboxButtonGroup(labels=["HIGH", "MODERATE", "LOW", "MODIFIER", "unknown"],
                                                active=[0, 1, 2, 3, 4])

    p_CustomJSFilter = functools.partial(CustomJSFilter, code="""
        const ld_start = slider.value[0]
        const ld_end = slider.value[1]
        const active_consq = button.active
        const indices = []

        for (var i=0; i < d.get_length(); i++) {
            var ld = d.data["ld"][i]
            var consq = d.data["ensembl_consequence_level"][i]
            if (ld_start <= ld && ld <= ld_end && active_consq.includes(consq)) {
                indices.push(i)
            }
        }
        return indices
    """)

    p_genome_circle = functools.partial(genome.circle, 
        'ps', 
        'logp', 
        line_width = 2, 
        size = 9, 
        fill_color = 'col', 
        line_color = "black", 
        line_alpha = 'col_assoc',
    )


    no_assoc_source = ColumnDataSource(e.loc[e['col_assoc']==0])
    filter_args = dict(slider = range_slider, button = checkbox_button, d = no_assoc_source)
    no_assoc_js_filter = p_CustomJSFilter(args = filter_args)
    no_assoc_view = CDSView(source = no_assoc_source, filters = [no_assoc_js_filter])
    p_genome_circle(source = no_assoc_source, view = no_assoc_view, legend_label = 'no assoc')

    assoc_source = ColumnDataSource(e.loc[e['col_assoc']==1])
    filter_args = dict(slider = range_slider, button = checkbox_button, d = assoc_source)
    assoc_js_filter = p_CustomJSFilter(args = filter_args)
    assoc_view = CDSView(source = assoc_source, filters = [assoc_js_filter])
    p_genome_circle(source = assoc_source, view = assoc_view, legend_label = 'assoc')
    
    # Make LD range slider
    range_slider.js_on_change('value', CustomJS(args=dict(no_assoc=no_assoc_source, assoc=assoc_source), code="""
    no_assoc.change.emit(); assoc.change.emit()
    """))
    checkbox_button.js_on_click(CustomJS(args=dict(no_assoc=no_assoc_source, assoc=assoc_source, ), code="""
    no_assoc.change.emit(); assoc.change.emit()
    """))

    genome.legend.location = "top_left"
    genome.legend.click_policy="hide"

    # Slightly extend the edges for viewability
    chrom = set(e['chrom']).pop()
    x_start = e['ps'].min()
    x_end = e['ps'].max()
    extend = (x_end - x_start) * 0.025
    genome.x_range.start = x_start - extend
    genome.x_range.end = x_end + extend

    y_start = e['logp'].min()
    y_end = e['logp'].max()
    extend = (y_end - y_start) * 0.025
    genome.y_range.start = y_start
    genome.y_range.end = y_end + extend

    genome.add_layout(Title(text="logp", align="center"), "left")

    ## Make GeneView plot
    geneview = GeneView(height = 500)
    geneview.update_view(genes)
    geneview.add_layout(Title(text="base position", align="center"), "below")
    geneview.x_range = genome.x_range

    # Add some metadata at the bottom of the plot
    geneview.add_layout(Title(text=f'chr{chrom}:{x_start}-{x_end}', align="right"), "below")

    timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')
    timestamp_text = f'Plot generated: {timestamp}'
    geneview.add_layout(Title(text=timestamp_text, align="right"), "below")
    geneview.add_layout(Title(text=f'Version: v{__version__}', align="right"), "below")

    # Add centromere region info on geneview plot
    cen_start, cen_end = _interactive_manh.get_centromere_region(chrom, build)
    cen_coordinate = set(range(cen_start, cen_end))
    region = set(range(x_start, x_end))
    cen_overlap = region.intersection(cen_coordinate)

    if (len(cen_overlap)>0): # Add indication of the centromeric region in the plot
        perc_overlap=int((len(cen_overlap)/len(region))*100)
        logger.info("\t\t\t    {0}% of the genomic region overlaps a centromere".format(perc_overlap))

        xs = min(cen_overlap)
        xe = max(cen_overlap)+1
        cen_median = max(cen_overlap)-int((max(cen_overlap)-min(cen_overlap))/2)
        geneview.segment(x0=xs, x1=xe, y0=0.5, y1=0.5, line_width=100, color="grey")
        cen = dict(x=[cen_median], y=[0.9], text=["Centromeric_region"])
        labels = LabelSet(x="x", y="y", text="text", source=ColumnDataSource(cen))
        geneview.add_layout(labels)


    ## Make peakplot
    peakplot = gridplot([[range_slider], [checkbox_button], [genome], [geneview]])
    return peakplot


def make_peakplot(infile, chrcol, pscol, a1col, a2col, pvalcol, mafcol, build, logger, vep_ld=0.1):
    ## Prepare data to create peakplot
    e, genes = make_view_data(infile, chrcol, pscol, a1col, a2col, pvalcol, mafcol, build, logger, vep_ld)
    peakplot = _create_peakplot(e, genes, build, logger)

    return peakplot


def output_peakplot(infile, outfile, title, pvalcol, pscol, mafcol, chrcol, a1col, a2col, build, logger, vep_ld=0.1):
    html = f'{outfile}.html'
    csv = f'{outfile}.csv'
    output_file(filename = html, title = title)

    e, genes = make_view_data(infile, chrcol, pscol, a1col, a2col, pvalcol, mafcol, build, logger, vep_ld)
    peakplot = _create_peakplot(e, genes, build, logger)
    
    e.to_csv(csv, header = True, index = False)
    save(peakplot)

