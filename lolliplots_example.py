# Lolliplots Example File
# -----------------------
lolliplot_folder = '/slade/home/user/scripts'
import sys
sys.path.append(lolliplot_folder)
import lolliplots as lp
import pandas as pd

lp.get_exon_data() # Run to download the ensemble exon data necessary for plots

## Important necessary files:
exon_file = "ensembl_exon_positions/b38_downloaded" #exon file downlaoded form get_exon_data
associations = "/slade/projects/UKBB/DNA_Nexus/bw_raw_raw_regenie_burden_dnanexus_2022-02-17/Single_Variant_bw_raw_Step2_Chr15_bw_raw.regenie" #Regenie variants output
masks = "/slade/projects/UKBB/DNA_Nexus/set_lists_450k_v2/annotations_chr15.txt" # Annotations file for gene
gene_name='IGF1R'
transcript_id='ENST00000649865'

exons = lp.read_exon_locs(exon_file) # Read exome locations in

selected_exons_df = exons.loc[exons['Transcript stable ID'] == transcript_id] # filter results to transcript
assoc_df = pd.read_csv(associations, delim_whitespace=True, comment='#') # Read regenie output
masks_df = pd.read_csv(masks, delim_whitespace=True, names=['ID','TID','MASK']) # Read masks file

# Merge gene results with masks (so that only selected exomes are used)
assoc_results = pd.merge(masks_df.loc[masks_df['TID'].str.contains(transcript_id),:], assoc_df, how='inner', on='ID' )

# Create lolliplot figure
fig1 = lp.lolliplot_raw(assoc_results, selected_exons_df, gene_name)

# Reduce gaps before creating figure (Note: two outputs from lp.reduce_gaps(), new asociation results file and new exons file)
assoc_short_gap, selected_exons_short_gap = lp.reduce_gaps(assoc_results, selected_exons_df)
fig2 = lp.lolliplot_raw(assoc_short_gap, selected_exons_short_gap, gene_name)

print('Example figure of exomes in "real" size')
fig1.show()
print('Example figure of exomes only (introns removed)')
fig2.show()

# Save the figure
# fig2.write_html('lolliplots_example.html')
