#!/usr/bin/env python
# coding: utf-8

# In[1]:


## Code based off Robin Beaumont's previous work on the Lollipop Plots that were created in Cambridge

import os as _os
import os.path as _ospath
import numpy as _np
import pandas as _pd
import subprocess as _subprocess
import plotly.express as _px
import plotly.graph_objects as _go


def get_exon_data():
    XML='''<?xml version="1.0" encoding="UTF-8"?>
	<!DOCTYPE Query>
	<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
				
		<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
			<Attribute name = "ensembl_gene_id" />
			<Attribute name = "ensembl_transcript_id" />
			<Attribute name = "exon_chrom_start" />
			<Attribute name = "exon_chrom_end" />
			<Attribute name = "start_position" />
			<Attribute name = "end_position" />
			<Attribute name = "external_gene_name" />
			<Attribute name = "ensembl_exon_id" />
		</Dataset>
	</Query>'''
    XML = XML.replace('\n','').replace('\t','')

    download_url = 'http://www.ensembl.org/biomart/martservice?query=' + XML
    file = 'ensembl_exon_positions/b38_downloaded'
    exists = _ospath.isfile('ensembl_exon_positions/b38_downloaded')

    if not exists:
        print('Ensembl Exon positions does not exist')
        repeat = True
        while repeat:
            print(f'Current folder: {_os.getcwd()}')
            download_yes = input('''Would you like to download them to folder: 
            ensemble_exon_positions/b38_downloaded
            (if folder doesn't exist, it will be created) [y/n]''')
            if download_yes.lower() == 'y':
                if not _ospath.exists('ensembl_exon_positions'):
                    _os.mkdir('ensembl_exon_positions')
                _subprocess.call(['wget','-O',file, download_url])
                repeat = False
            elif download_yes.lower() == 'n':
                repeat = False
            else:
                print('Input not understood please input "y" or "n" only')

def read_exon_locs(file):
    HEADERS = ['Gene stable ID',
           'Transcript stable ID',
           'Exon start (bp)',
           'Exon stop (bp)',
           'Gene start (bp)',
           'Gene end (bp)',
           'Gene name',
           'Exon Name']
    exons = _pd.read_csv(file, delim_whitespace=True, names=HEADERS)
    
    return exons

def lolliplot_raw(results, exon_info,
                  title,
                  fig_height = 500,
                  fig_width = 1400,
                  ex_start_col = 'Exon start (bp)',
                  ex_stop_col = 'Exon stop (bp)',
                  lolli_x = 'GENPOS',
                  lolli_y = 'LOG10P',
                  lolli_size = 'BETA',
                  lolli_col = 'MASK',
                  lolli_direction = 'BETA',
                  lollipop_max_size = 5.,
                  lollipop_stem_width = 0.1):
    '''
    Makes a lolliplot based on the given results file and exon lcoations file (unmanipulated)
    
    Parameters
    ----------
    results: pd.DataFrame
    
    exon_info: pd.DataFrame
    
    title: Str
    
    Defaults
    --------
    fig_height: int: 500
        Figure height (can be adjusted using fig.update_layout())
    fig_width: int: 1400
        Figure height (can be adjusted using fig.update_layout())
    ex_start_col: Str: 'Exon start (bp)'
        Name of column of exon_info of exome start locations
    ex_stop_col: Str: 'Exon stop (bp)'
        Name of column of exon_info of exome stop locations
    lolli_x: Str: 'GENPOS'
        Name of column of results df for x axis
    lolli_y: Str: 'LOG10P'
        Name of column of results df for y axis
    lolli_size: Str: 'BETA'
        Name of column of results df for size of bubbles
    lolli_col: Str: 'MASK'
        Name of column of results df for colour of bubbles
    lolli_direction: Str: 'BETA'
        Name of column of reults df for direction of bubbles
    lollipop_max_size: float: 5.
        Maximum size of bubbles
    lollipop_stem_width: float: 0.1
        Width of the lollipop stems
    
    Returns
    -------
    Plotly figure
    '''
    
    # Make the Exome Rectangle Points
    exon_x = []
    exon_y = []
    for row in exon_info.index:
        exon_x += 2*[exon_info[ex_start_col][row]]
        exon_x += 2*[exon_info[ex_stop_col][row]]
        exon_x += [exon_info[ex_start_col][row]]
        exon_x += [None]
        exon_y += [0.1, -0.1, -0.1, 0.1, 0.1, None]
    
    # Draw the figure and rectangles
    fig = _go.Figure(
        _go.Scatter(x=exon_x,
                   y=exon_y,
                   mode='lines',
                   fill='toself',
                   showlegend=False,
                   name='Exon regions',
                   line=dict(color='black', width=0.5),
                   fillcolor='rgba(100,100,255,1)')
    )
    
    # Draw lollipops onto the exomes 
    for color in results[lolli_col].unique():
        df_filt = results.loc[results[lolli_col]==color, :]
        
        # Create data to draw the lines to the exome
        lines_x = []
        lines_y = []
        for row in df_filt.index:
            lines_y += [df_filt[lolli_y][row] * _np.sign(df_filt[lolli_direction][row]), 0, None]
            lines_x += [df_filt[lolli_x][row], df_filt[lolli_x][row], None]

        # Add bubble scatter and lines
        fig.add_trace(
            _go.Scatter(
                x = df_filt[lolli_x], 
                y = df_filt[lolli_y]*_np.sign(df_filt[lolli_direction]),
                mode='markers',
                marker=dict(size=abs(df_filt[lolli_size]),
                            sizeref=2.*results[lolli_size].max()/(lollipop_max_size**2)),
                name=color,
                legendgroup=color)
        )
        
        fig.add_trace(
            _go.Scatter(
                x = lines_x,
                y = lines_y,
                mode='lines',
                showlegend=False,
                line=dict(color='black', width=lollipop_stem_width),
                hoverinfo='skip',
                name = color,
                legendgroup = color)
        )
        
    # Update the figure so legend colour is visible
    fig.update_layout(legend= {'itemsizing': 'constant'})
    
    # Reverse the order so that exome locations and the bubbles are on top
    fig.data = fig.data[::-1]
    
    # Add titles
    fig.update_layout(
        title=title,
        xaxis_title='Exon Positions',
        yaxis_title=lolli_y,
        legend_title=lolli_col,
        width=fig_width,
        height=fig_height
    )
    return fig

def reduce_gaps(gene_results, exon_info, new_gap=10,
                ex_start_col = 'Exon start (bp)',
                ex_stop_col = 'Exon stop (bp)',
               gene_pos_col = 'GENPOS'):
    '''
    Updated the results and exon dataframes to reduce the gap between exons
    
    Parameters
    ----------
    gene_results: pd.DataFrame
    
    exon_info: pd.DataFrame
    
    Defaults
    --------
    new_gap: int: 10

    ex_start_col: Str: 'Exon start (bp)'
        Name of column of exon_info of exome start locations
    ex_stop_col: Str: 'Exon stop (bp)'
        Name of column of exon_info of exome stop locations
    geen_pos_col: Str: 'GENPOS'
        Name of colmumn of geen positions
    
    Returns
    -------
    Tuple:
    1. Updated gene_results dataframe
    2. Updated exon_info dataframe
    '''
    # Calculate gap distances
    exon_sorted = exon_info.sort_values(ex_start_col).reset_index(drop=True)
    gap = _np.roll(exon_sorted[ex_start_col], shift=-1) - exon_sorted[ex_stop_col]
    gap = gap[:-1]
    
    results_out = gene_results.copy()
    both_col = [ex_start_col, ex_stop_col]
    for i in range(len(gap)):
        exon_sorted.loc[i+1:, both_col] = exon_sorted.loc[i+1:, both_col].applymap(lambda x: x - gap[i] + new_gap)
        my_min = exon_sorted.loc[i+1, ex_start_col]
        results_out.loc[(results_out[gene_pos_col] > my_min), gene_pos_col] = results_out.loc[(results_out[gene_pos_col] > my_min), gene_pos_col].map(lambda x: x - gap[i] + new_gap)
        
    return results_out, exon_sorted

def example():
    ## Important necessary files:
    exon_file = "ensembl_exon_positions/b38_downloaded"
    associations = "/slade/projects/UKBB/DNA_Nexus/bw_raw_raw_regenie_burden_dnanexus_2022-02-17/Single_Variant_bw_raw_Step2_Chr15_bw_raw.regenie"
    masks = "/slade/projects/UKBB/DNA_Nexus/set_lists_450k_v2/annotations_chr15.txt"
    gene_name='IGF1R'
    transcript_id='ENST00000649865'
    #max_af=0.001
    
    exons = read_exon_locs(exon_file)
    
    selected_exons_df = exons.loc[exons['Transcript stable ID'] == transcript_id]
    assoc_df = _pd.read_csv(associations, delim_whitespace=True, comment='#')
    masks_df = _pd.read_csv(masks, delim_whitespace=True, names=['ID','TID','MASK'])
    
    # Merge gene results with masks (so that only selected exomes are used)
    assoc_results = _pd.merge(masks_df.loc[masks_df['TID'].str.contains(transcript_id),:], assoc_df, how='inner', on='ID' )
    
    fig1 = lolliplot_raw(assoc_results, selected_exons_df, gene_name)
    assoc_short_gap, selected_exons_short_gap = reduce_gaps(assoc_results, selected_exons_df)
    fig2 = lolliplot_raw(assoc_short_gap, selected_exons_short_gap, gene_name)
    
    print('Example figure of exomes in "real" size')
    fig1.show()
    print('Example figure of exomes only (introns removed)')
    fig2.show()

if __name__ == "__main__":
    example()