
# import required packages
import pandas as pd
import os
import sys
import clumpsptm
# import pymol
from tqdm import tqdm
import glob
from typing import Union
import matplotlib.pyplot as plt
import numpy as np
import ast

# import general PANOPLY utilities
import argparse
import yaml

# import utility functions
import tarfile
import warnings
import random
import re
import pprint
from datetime import datetime

# # import debugging functions
# import inspect # for debugging
# import shutil # for file copy



##################################
####   Parameter Management   ####
##################################

# create custom illegal-argument error, for use when checking parameters
class IllegalArgumentError(ValueError):
    pass

# # configure ProDy to not printout messages, to avoid
# from prody import confProDy
# confProDy(verbosity='error')

### import command line parameters
parser = argparse.ArgumentParser(prog = 'ClumpsPTM Postprocessing', description="Script for processing the results of ClumpsPTM")

parser.add_argument("-r", "--results_tar", type=str, help="Directory with results from ClumpsPTM analysis", required=True)
# parser.add_argument("-m", "--mapping_file", type=str, help="ClumpsPTM mapping file with mappings to PDB archive", required=True)

parser.add_argument("-f", "--fdr_threshold", type=float, help="Threshold for minimum theoretical p-value to consider for FDR")

# parser.add_argument("-i", "--accession_col", type=str, help="GCT rdesc column with accession IDs. Must match the ID type of the provided FASTA file.")
# parser.add_argument("-v", "--variable_sites_col", type=str, help="GCT rdesc column with PTM variable site(s) (e.g. 'T527t')")
# parser.add_argument("-s", "--variable_sites_sep", type=str, help="Separator for variable sites (e.g. ' ' is the separator for 'T972t S977s')")
# parser.add_argument("-g", "--gene_column", type=str, help="GCT rdesc column with HUGO Gene Symbols")


parser.add_argument("-o", "--output_prefix", type=str, help="Output prefix to prepend to output filenames.", default="results")
parser.add_argument("-y", "--yaml", type=str, help="Path to .yaml file with parameters.", required=True)


# import from command line
args = parser.parse_args()

# # optional testing args
# args = parser.parse_args(
#     ["-r" "/opt/input/ODG_v3_NMF.consensus.core.k3_clumps_runs.tar", \
#      "-f" "0.1", \
#      # "-i" "id.description", \
#      # "-s" " ", \
#      # "-v" "variableSites", \
#      # "-g" "geneSymbol", \
#      "-o" "ODG_v3_NMF.consensus.core.k3", \
#      "-y" "/opt/input/master-parameters.yaml"]
# )


####   Default Parameter Import   ####

with open(args.yaml, 'r') as file:
    yaml_dict = yaml.safe_load(file)

# # override missing parameters with yaml defaults
# if (args.accession_col==None):
#     # args.accession_col = yaml_dict['global_parameters']['gene_mapping']['protein_id_col']
#     args.accession_col = yaml_dict['panoply_ptm_normalization']['accession_number_colname']

# if (args.variable_sites_col==None):
#     args.variable_sites_col = yaml_dict['panoply_clumps_ptm']['mapping']['variable_sites_col']

# if (args.variable_sites_sep==None):
#     args.variable_sites_sep = yaml_dict['panoply_clumps_ptm']['mapping']['variable_sites_sep']

# if (args.gene_column==None):
#     args.gene_column = yaml_dict['global_parameters']['gene_mapping']['gene_id_col']



# print parameters
print('\n\nParameters:')
pprint.pp(args.__dict__)
print('\n')



# tmp params
# ome_types = ['acetylome','ubiquitylome','phosphoproteome']
# ome_types = ['acetylome','phosphoproteome']


####################################
####   File & Directory Setup   ####
####################################

# set up output directories
os.makedirs(args.output_prefix, exist_ok=True)

out_dir_figs=os.path.join(args.output_prefix,'figures')
os.makedirs(out_dir_figs, exist_ok=True)
out_dir_pymol = os.path.join(out_dir_figs,'pymol')
os.makedirs(out_dir_pymol, exist_ok=True)
out_dir_dotplots = os.path.join(out_dir_figs,'dotplots')
os.makedirs(out_dir_dotplots, exist_ok=True)



# extract tarfile contents
tar = tarfile.open(args.results_tar)

tarfile_outdir = args.output_prefix
tar.extractall(tarfile_outdir)


#####################################
####   Data Import / Wrangling   ####
#####################################



results_list = list()

results_dirs = [os.path.basename(fn) for fn in glob.glob(os.path.join(tarfile_outdir,"**"))]

# import TSV results files
for dir in tqdm(results_dirs):
    for file in glob.glob(os.path.join(tarfile_outdir, dir, "*.tsv")):
        _df = pd.read_csv(file, sep='\t', index_col=0)
        _df['id'] = dir # use directory as analysis ID
        _df['subval'] = re.sub('^(.+)_(.+?)_results$', '\\1', dir) # use directory as analysis ID
        _df['direction'] = re.sub('^(.+)_(.+?)_results$', '\\2', dir) # use directory as analysis ID
        results_list.append(_df)

results_df = pd.concat(results_list)



# add FDR correction
res = list()

for idx in tqdm(np.unique(results_df['id'])):
    _df = results_df[results_df['id']==idx].drop(columns=['fdr_max_pval','fdr_pass','fdr_corr'])
    _df = clumpsptm.utils.add_corrected_fdr(_df, args.fdr_threshold) # has a noisy printout
    res.append(_df)

results_df = pd.concat(res)



res = list()

for idx in tqdm(np.unique(results_df['id'])):
    _df = results_df[results_df['id']==idx].drop(columns=['fdr_max_pval','fdr_pass','fdr_corr'])
    _df = clumpsptm.utils.add_corrected_fdr(_df, args.fdr_threshold, weight_thresh_by_n=True) # has a noisy printout
    res.append(_df)

#results_01_weight_df
results_fdrThresh_weight_df = pd.concat(res)



# TODO: the variableSite column gets exploded in the original pass-- and we don't save the original values. Each PTM will appear as if its a single-site PTM in the mapping file.

# # Highlight any hits that are derived from a singular peptide
# df = pd.read_csv(args.mapping_file, sep='\t', index_col=0)
# acc_var_df = df[[args.accession_col,args.variable_sites_col]].reset_index().set_index(args.accession_col).drop_duplicates()
# var_site_col = 'tmp'
# acc_var_df[var_site_col] = acc_var_df[args.accession_col].map(nchar)
# acc_var_df[variable_sites_col] = [str.split(args.variable_sites_sep) for str in acc_var_df[args.variable_sites_col]] # convert variableSites column into list


# def traceback_peptide3s(row):
#     """Traceback peptide 3s"""
    
#     if row['clumpsptm_input_n'] != 3:
#         return None
#     else:
#         vs3 = {x if len(x.split(ptm_split)) == 3 else None for x in row['variableSites']}.pop()
#         print(vs3)
#         try:
#             return acc_var_df[acc_var_df['variableSites']==vs3].loc[row.name]['id']
#         except:
#             return None

# results_df['trace3s'] = results_df.apply(traceback_peptide3s,1).values





###############################
####   Plotting Dotplots   ####
###############################



def plot_pair(group, results_df, n_to_plot=20):
    """Plot paring."""
    ome_types = np.unique(results_df['clumpsptm_sampler']) # get unique feature-types
    if ('ptm' in ome_types): # if we ran combined
        np.insert(np.delete(ome_types, ome_types=='ptm'), 0, 'ptm') # move 'ptm' to beginning

    fig, axes = plt.subplots(2, len(ome_types), figsize=(10,14)) # TODO: needs to be adjusted if we have additional omes
    # for each feature
    for j,feature in enumerate(ome_types):
    	# plot positive/negative features
        for i,direction in enumerate(['negative','positive']):
            _df = results_df[results_df['id']=="{}_{}_results".format(group, direction)]
            clumpsptm.vis.dotplot(
                _df.loc[_df[_df['clumpsptm_sampler']==feature].sort_values(
                    "clumpsptm_pval").iloc[:n_to_plot,:].index], # NOTE: the final .index grabs ALL features associated with the top ID, so you can a feature in the context of the other omes 
                sort_by=feature,
                x='clumpsptm_pval',
                ax=axes[i,j],
                thresh=args.fdr_threshold,
            )
            # add legend
            axes[i,j].legend().remove()
            axes[i,j].set_yticklabels(axes[i,j].get_yticklabels(), fontsize=14)
            axes[i,j].set_xlabel(r"$-log_{10}$ p-value", fontsize=16)
            # add positive/negative labels
            if direction=='positive': direction = "(+)"
            if direction=='negative': direction = "(-)"
            # set title and y-ticks
            axes[i,j].set_title("{} {}".format(feature.capitalize(), direction), fontsize=18)
            axes[i,j].set_yticklabels(axes[i,j].get_yticklabels(), fontsize=14)
    # plot subplots
    plt.suptitle(group, x=0.55, y=1.025, fontsize=20)
    plt.tight_layout()




for group in np.unique(results_df['subval']):
    plot_pair(group, results_df)
    plt.savefig(os.path.join(out_dir_dotplots, "{}_dotplot.pdf".format(group)), dpi=300, bbox_inches='tight')




####################################
####   Plotting PyMol Figures   ####
####################################


## Creation of Pymol Files
## --------------------------

# for group in np.unique(results_df['id']):
#     os.makedirs(os.path.join('figures/pymol', group), exist_ok=True)
    
#     for feat in np.unique(results_df['clumpsptm_sampler']):
#         _df = results_df[(results_df['id']==group) & 
#                          (results_df['clumpsptm_sampler']==feat)
#                         ].sort_values('clumpsptm_pval').reset_index()
#         _df.index = _df.index.astype(str)
#         _out_dir = os.path.join('figures/pymol',group,feat)
#         clumpsptm.vis.create_pymols_from_result(_df, out_dir=_out_dir, include_idx_in_name=True)







#############################
####   Sumamry Figures   ####
#############################

# TODO: decide if we should subset this or not
to_plot = np.unique(results_df['id'])


# ### Summary Figure
results_filt_df = results_df[results_df['id'].isin(to_plot)]

counts_df = results_filt_df.groupby(['id','clumpsptm_sampler']).size().reset_index().set_index(["id","clumpsptm_sampler"])
counts_df.columns = ['NS']

counts_pval_df = results_df[results_df['clumpsptm_pval']<args.fdr_threshold].groupby(['id','clumpsptm_sampler']).size().reset_index().set_index(["id","clumpsptm_sampler"])
counts_pval_df.columns = ['< {} P-Value'.format(args.fdr_threshold)]

counts_fdr_df = results_df[results_df['fdr_corr']<args.fdr_threshold].groupby(['id','clumpsptm_sampler']).size().reset_index().set_index(["id","clumpsptm_sampler"])
counts_fdr_df.columns = ['< {} FDR'.format(args.fdr_threshold)]

counts_df = counts_df.join(counts_pval_df).join(counts_fdr_df).fillna(0).astype(int).sort_values('NS')
counts_df = counts_df.reset_index()
#counts_df = counts_df[counts_df['id'].str.contains('positive') | counts_df['id'].str.contains('negative')].set_index(['id','clumpsptm_sampler'])

# pick -ome to sort data by
ome_types = np.unique(counts_df['clumpsptm_sampler'])
if ('ptm' in ome_types):
    np.insert(np.delete(ome_types, ome_types=='ptm'), 0, 'ptm') # move 'ptm' to beginning

# get _order using first PTM type
_counts_df = counts_df.reset_index()
_counts_df = _counts_df[_counts_df['clumpsptm_sampler']==ome_types[0]].set_index("id") # TODO: use first PTM for now, but eventually use combined dataset
# _order = _counts_df.sort_values(by=['NS']).index 
_order = _counts_df.sort_values(by=['id']).index # sort in order for now


### Summary Figure
## Create figure

fig,axes = plt.subplots(1,len(ome_types),figsize=(12,5),sharey=False)

for i,ome in enumerate(ome_types):
	_counts_df = counts_df.reset_index()
	_counts_df = _counts_df[_counts_df['clumpsptm_sampler']==ome].set_index("id").loc[_order][['NS', '< {} P-Value'.format(args.fdr_threshold), '< {} FDR'.format(args.fdr_threshold)]]
	_counts_df.plot(kind='barh', stacked=True, ax=axes[i], linewidth=1, width=0.8, edgecolor='black', color=['lightgrey','orange','red'])
	# set Axis labels
	axes[i].set_title(ome, fontsize=16)
	axes[i].set_ylabel("")
	axes[i].set_xlabel("# Proteins", fontsize=16)
	axes[i].legend().remove()
	if (i>0): 
		axes[i].set_yticks([])

# add final legend
axes[i].legend(loc='center left', bbox_to_anchor=(1, 0.5))

axes[0].set_yticklabels([x.get_text().replace("_", " ").upper().replace("POS RESULTS","(+)").replace("NEG RESULTS","(-)")  for x in axes[0].get_yticklabels()])
axes[0].set_xlim([0,1200])
axes[1].set_xlim([0,1200])
axes[2].set_xlim([0,1200])


plt.tight_layout()
plt.savefig(os.path.join(out_dir_figs,"clumpsptm_summary_figure.pdf"), dpi=300, bbox_inches='tight')



