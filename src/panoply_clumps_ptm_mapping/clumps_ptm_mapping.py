#### Mapping Script for ClumpsPTM

import clumpsptm
import pandas as pd
import os
import numpy as np
import glob
from tqdm import tqdm
import subprocess
import ast
from agutil.parallel import parallelize2
#import matplotlib.pyplot as plt

#import requests, sys # for FASTA import
import re # for regular expressions
import itertools
import argparse
import yaml


# import GCT management tools
import cmapPy
from cmapPy.pandasGEXpress.parse_gct import parse
from cmapPy.pandasGEXpress.write_gct import write

# import utility functions
import warnings
import random
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

# create str2bool() for importing args.run_combined parameter
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# configure ProDy to not printout messages, to avoid
from prody import confProDy
confProDy(verbosity='error')

### import command line parameters

parser = argparse.ArgumentParser(prog = 'ClumpsPTM Mapping', description="Script for mapping to PDBs for ClumpsPTM")
parser.add_argument("-p", "--phosphoproteome_gct", type=str, help="Phosphoproteome GCT file")
parser.add_argument("-a", "--acetylome_gct", type=str, help="Acetylome GCT file")
parser.add_argument("-u", "--ubiquitylome_gct", type=str, help="Ubiquitylome GCT file")

parser.add_argument("-f", "--FASTA_ref_file", type=str, help="Reference FASTA file with all relevant sequences for your dataset. Used to BLAST against Uniprot sequences.", required=True)
parser.add_argument("-t", "--FASTA_sep_type", type=str, help="Separator for FASTA sequences in FASTA reference file. Supported values are \"cptac\" for \" \", or \"gencode\" for \"|\".')")

parser.add_argument("-i", "--accession_col", type=str, help="GCT rdesc column with accession IDs. Must match the ID type of the provided FASTA file.")
parser.add_argument("-g", "--gene_column", type=str, help="GCT rdesc column with HUGO Gene Symbols")

# PTM Site Management
parser.add_argument("-v", "--variable_sites_col", type=str, help="GCT rdesc column with PTM variable site(s) (e.g. 'T527t')")
parser.add_argument("-s", "--variable_sites_sep", type=str, help="Separator for variable sites (e.g. ' ' is the separator for 'T972t S977s')")
parser.add_argument("--keep_multi_sites", type=str2bool, help="Should multi-site PTMs be mapped?.")
parser.add_argument("--filter_duplicate_sites", type=str2bool, help="Should multi-site PTMs be filtered to remove sites that were observed as single-sites?")

parser.add_argument("-b", "--PDB_DIR", type=str, help="Directory with PDB structures.", required=True)
parser.add_argument("--UNIPROT_SWISSPROT", type=str, help="Reference FASTA file with all relevant UNIPROT sequences, to BLAST your sequences to.", required=True)
parser.add_argument("--SIFTS_DB", type=str, help="SIFTS database containing mapping between UNIPROT IDs and PDB IDs.", required=True)

parser.add_argument("-o", "--output_prefix", type=str, help="Output prefix to prepend to output filenames.", default="")
parser.add_argument("-y", "--yaml", type=str, help="Path to .yaml file with parameters.", required=True)

parser.add_argument("-n", "--num_threads", type=int, help="Number of threads available", default = 1) # assume 1 thread if not provided, just to be extra safe
parser.add_argument("-d", "--DEBUG_MODE", help="Run in debugging mode; limit the number of FASTA files processed.", action = "store_true")
parser.add_argument("--DEBUG_RNG", type=int, help="RNG seed for debugging.")


# import from command line
args = parser.parse_args() # import from command line

# # testing arguments manually
# args = parser.parse_args(["-p" "/opt/input/phosphoproteome-subset.gct", \
#     "-a" "/opt/input/acetylome-subset.gct", \
#     "-u" "/opt/input/ubiquitylome-subset.gct", \
#     "-f" "/opt/input/Gencode_v39.basicPCnr2.642contams.fasta", \
#     "--PDB_DIR", "/pdbs", \
#     "--UNIPROT_SWISSPROT", "/opt/input/uniprot_sprot.fasta", \
#     "--SIFTS_DB", "/opt/input/pdb_chain_uniprot.tsv", \
#     "-i", "id.description", \
#     "-o", "ODG_v3", \
#     "-y", "/opt/input/master-parameters.yaml",\
#     "--DEBUG_MODE", \
#     "-n", "12"])
# # args = parser.parse_args("--phosphoproteome_gct /opt/input/var_map_fromOutput_phosphoproteome.gct --acetylome_gct /opt/input/var_map_fromOutput_acetylome.gct --PDB_DIR /pdbs --FASTA_ref_file /opt/input/RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams_553smORFs.fasta --FASTA_sep_type cptac --accession_col accession_number --variable_sites_col variableSites_edited --UNIPROT_SWISSPROT /opt/input/Freeze_061721_clumpsptm_ref_uniprot_uniprot_sprot.fasta --SIFTS_DB /opt/input/Freeze_061721_clumpsptm_ref_uniprot_pdb_chain_uniprot.tsv --output_prefix pancan --yaml /opt/input/master-parameters.yaml --num_threads 8 --DEBUG_MODE".split())
# # args = parser.parse_args("-p /opt/input/ODG-v3-phosphoproteome-SpectrumMill-ratio-QCfilter-NArm.gct -u /opt/input/ODG-v3-ubiquitylome-SpectrumMill-ratio-QCfilter-NArm.gct -a /opt/input/ODG-v3-acetylome-SpectrumMill-ratio-QCfilter-NArm.gct --PDB_DIR /pdbs -f /opt/input/Ensembl.human.hg19.clean3nr.602contams_20230913.fasta --UNIPROT_SWISSPROT /opt/input/uniprot_sprot.fasta --SIFTS_DB /opt/input/pdb_chain_uniprot.tsv -i id.description -o ODG_v3 -y /opt/input/master-parameters.yaml -n 12".split())



# ensure that at least one GCT has been provided
if (args.phosphoproteome_gct==None and \
    args.acetylome_gct==None and \
    args.ubiquitylome_gct==None):
    raise IllegalArgumentError("No GCT files provided. Please provide at least one PTM GCT file.")



### import default parameters from YAML

with open(args.yaml, 'r') as file:
    yaml_dict = yaml.safe_load(file)

# override missing parameters with yaml defaults

if (args.accession_col==None):
    # args.accession_col = yaml_dict['global_parameters']['gene_mapping']['protein_id_col']
    args.accession_col = yaml_dict['panoply_ptm_normalization']['accession_number_colname']

if (args.gene_column==None):
    args.gene_column = yaml_dict['global_parameters']['gene_mapping']['gene_id_col']

if (args.variable_sites_col==None):
    args.variable_sites_col = yaml_dict['panoply_clumps_ptm']['mapping']['variable_sites_col']

if (args.keep_multi_sites==None):
    args.keep_multi_sites = yaml_dict['panoply_clumps_ptm']['mapping']['keep_multi_sites']

if (args.filter_duplicate_sites==None):
    args.filter_duplicate_sites = yaml_dict['panoply_clumps_ptm']['mapping']['filter_duplicate_sites']

if (args.FASTA_sep_type==None):
    args.FASTA_sep_type = yaml_dict['panoply_clumps_ptm']['mapping']['FASTA_sep_type']

if (args.DEBUG_RNG==None):
    args.DEBUG_RNG = int(datetime.now().strftime("%H%M%S")) # set seed to current time

# ensure that FASTA_sep_type is a valid value
if (args.FASTA_sep_type!='gencode' and \
    args.FASTA_sep_type!='cptac'):
    raise IllegalArgumentError("FASTA_sep_type must be either 'gencode' (for '|') or 'cptac' (for ' '); the value '"+args.FASTA_sep_type+"' is not allowed")


# print parameters
print('\n\nParameters:')
pprint.pprint(args.__dict__)
print('\n')



####################################
####   File & Directory Setup   ####
####################################

# Downloaded PDB Database
# PDB_DIR = "/pdbs" # directory with PDB database
PDB_DIR = args.PDB_DIR # directory with PDB database

# Reference & Output Directories
REF_DIR = "/reference_files"
# REF_DIR = "/opt/input"
FASTA_DIR = "fasta_files"
OUT_DIR = "output_files"
os.makedirs(os.path.join(REF_DIR, FASTA_DIR), exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# save parameters file
with open(os.path.join(OUT_DIR,'params.yaml'), 'w') as f:
    # json.dump(args.__dict__, f)
    yaml.dump(args.__dict__, f) # dump to YAML


# GCT dictionary with all provided PTM GCTs
gcts_zipped = zip(['phosphoproteome', 'acetylome', 'ubiquitylome'], \
    [ args.phosphoproteome_gct, args.acetylome_gct, args.ubiquitylome_gct ])
gct_dict = {label: gct for label, gct in gcts_zipped if gct is not None}

# # Reference Files
# UNIPROT_SWISSPROT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
# UNIPROT_SWISSPROT = "uniprot_sprot.fasta"
UNIPROT_SWISSPROT = args.UNIPROT_SWISSPROT
# SIFTS_DB = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
# SIFTS_DB = "pdb_chain_uniprot.tsv"
SIFTS_DB = args.SIFTS_DB

# Other Parameters
n_threads = args.num_threads
accn_col = args.accession_col
ptm_col = args.variable_sites_col
ptm_split = args.variable_sites_sep
keep_multi_sites = args.keep_multi_sites
filter_duplicate_sites = args.filter_duplicate_sites
gene_col = args.gene_column
accn_type = 'inputID' # label to use for input ID (e.g. 'ENSEMBL' or 'refseq'). not really used; currently a placeholder.
unique_id_col = 'rid' # hardcode unique_id_col to "rid" since we're importing a GCT file. If we were using a pre-generated var_site.tsv file, this might be "id" instead.

DEBUG_MODE = args.DEBUG_MODE # running in DEBUG_MODE will limit the number of FASTA processed to DEBUG_N_PROT
DEBUG_N_PROT = 50

output_prefix = args.output_prefix+"_"



if (DEBUG_MODE): # if we're running DEBUG mode, warn the user
    DEBUG_RNG = np.random.RandomState(args.DEBUG_RNG) # create random state from seed
    print(f"#### WARNING: DEBUG_MODE is currently toggled ON. This setting should NOT be used for full-runs of Clumps-PTM")




################################################
####   Import GCTs / Create Feature File    ####
################################################

# import GCT(s)


rdesc_list = []

for label, gct_fn in gct_dict.items():
    try:
        # read in file
        _gct = parse(gct_fn)
        _rdesc = _gct.row_metadata_df # Access row and column metadata
        # add PTM type as a 'feature' column
        _rdesc['feature'] = label 
        # check for required columns
        try: 
            col=accn_col; _tmp=_rdesc[col]
            col=ptm_col; _tmp=_rdesc[col]
            col=gene_col; _tmp=_rdesc[col]
        except KeyError: # error if missing
            print(f"Missing rdesc column '{col}' in {label} GCT.".format(col=col,label=label))
        # append to rdesc list
        rdesc_list.append(_rdesc) # append to list
    except FileNotFoundError:
        print(f"Error: The file '{gct}' was not found.")


# combine relevant parts of rdesc
pmap_df = pd.concat(rdesc_list)

# write to Feature File for future use
FEATURE_FILENAME = os.path.join(OUT_DIR,output_prefix+"var_sites_combined.tsv")
pmap_df.to_csv(FEATURE_FILENAME, sep='\t')


###########################
####    File Import    ####
###########################
# # Import Reference Files
# clumpsptm.mp.dl_ref(REF_DIR, [UNIPROT_SWISSPROT, SIFTS_DB])

# Import PTM-Site Annotations & get list of accession numbers
pmap_df = pd.read_csv(FEATURE_FILENAME, sep='\t', index_col=0) # read back in feature_file. redundant but whatever.
#pmap_df['accn_base'] = [re.sub("\\.\d+?","",accn) for accn in pmap_df['id.description']] # add base accession number to pmap_df. only needed for downloading FASTAs.
accn_arr = pmap_df[accn_col].drop_duplicates() # get unique accession numbers

# Import Sifts Data-Base
sifts_df = pd.read_csv(SIFTS_DB, comment="#", sep='\t', low_memory=False)


print("#### FILE IMPORT COMPLETE  --- ", datetime.now())

##################################
####   Read in FASTA Files    ####
##################################

# split FASTAs into individual files
clumpsptm.mp.split_fastas(
    args.FASTA_ref_file,
    os.path.join(REF_DIR, FASTA_DIR),
    naming=args.FASTA_sep_type # determines which filename-splitting delimiter is used; "cptac"=" ", "gencode"="|"
)
individual_fastas_all = glob.glob(os.path.join(REF_DIR, FASTA_DIR, "*"))

# filter FASTAs to accession numbers in our dataset
#filt_in_df = [any(map(fasta.__contains__, accn_arr)) for fasta in individual_fastas_all]
pattern = re.compile('|'.join(map(re.escape, accn_arr)))
filt_in_df = [bool(pattern.search(fasta)) for fasta in individual_fastas_all] # slightly faster than map approach
#itertools.compress(individual_fastas_all,filt_in_df)
individual_fastas = list(itertools.compress(individual_fastas_all,filt_in_df))


print("\n#### FASTA FILE-SPLITTING COMPLETE  --- ", datetime.now())

if (DEBUG_MODE): # if we're running DEBUG mode, only use a subset of proteins
    print(f"#### WARNING: DEBUG_MODE is currently toggled ON. This setting should NOT be used for full-runs of Clumps-PTM")
    individual_fastas = DEBUG_RNG.choice(individual_fastas, DEBUG_N_PROT)



# # alternatively: import FASTA files based on accession numbers
# os.system('pip install requests')
# import requests
# import random
# server = "https://rest.ensembl.org"
# #ext = "/sequence/id/ENSG00000157764?type=protein;multiple_sequences=1"
# #ext = "/sequence/id/{accn}?type=protein;multiple_sequences=1".format(accn = "ENSP00000296755.7")
# #ext = "/sequence/id/{accn}?type=protein;format=fasta;multiple_sequences=1;species=human".format(accn = "ENSP00000296755")
# ext_pattern = "/sequence/id/{accn}?type=protein;format=fasta;multiple_sequences=1"

# n_proteins = 50
# for accn_full in accn_arr[random.sample(range(0,len(accn_arr)), n_proteins)]:
# 	accn_base = re.sub("\\.\d+?","",accn_full)
# 	ext = ext_pattern.format(accn=accn_base)
# 	r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"}) 
# 	if not r.ok:
# 		r.raise_for_status()
# 		sys.exit()
# 	entries = r.text.split(">")[1:] # split FAST into individual acceession numbers (pruning first, blank entry)
# 	for entry in entries:
# 		accn = entry.strip().split("\n")[0] # pull out accession number
# 		fn = "{dir}/{subdir}/{accn}.seq".format(accn=accn, dir=REF_DIR, subdir=FASTA_DIR)
# 		with open(fn, "w") as f:
# 			f.write(r.text)

# individual_fastas = glob.glob(os.path.join(REF_DIR, FASTA_DIR, "*"))
# #individual_fastas[:10]
# if len(individual_fastas) == 0:
# 	raise Exception("No FASTA files available for BLAST")
# os.listdir(os.path.join(REF_DIR, FASTA_DIR))






##################################
####    BLAST FASTA Files     ####
##################################
BLAST_DIR = accn_type+"_to_uniprot_blast"
clumpsptm.mp.blast_sequences(
    UNIPROT_SWISSPROT,
    individual_fastas,
    output_dir=REF_DIR,
    n_threads=n_threads,
    db_title="UniprotDB",
    blast_dir_name=BLAST_DIR,
    seq_db_name="uniprot_db"
)
print("#### BLAST COMPLETE  --- ", datetime.now())

# collect BLASTed files
blasted_files_all = glob.glob(os.path.join(REF_DIR, BLAST_DIR, "*")) # find all files in BLAST directory
# check for empty blasted_files so clumpsptm.mp.get_blast_hits_with_sifts() doesn't fail
blasted_files = [f for f in blasted_files_all if os.path.getsize(f)!=0]


#blasted_files[:10]
if len(blasted_files) == 0:
	raise Exception("No BLASTed files...")

# import BLAST results that have SIFTS hits
#mapped_acc_df = clumpsptm.mp.get_blast_hits_with_sifts(blasted_files, sifts_df)
mapped_acc_df = clumpsptm.mp.get_blast_hits_with_sifts(blasted_files, sifts_df, filter_only_human=True)
#mapped_acc_df
if mapped_acc_df.shape[0]==0:
	raise Exception("No accession numbers with SIFTS hits")


# filter sifts to those in the BLAST results
sifts_filt_df = sifts_df[sifts_df['SP_PRIMARY'].isin(mapped_acc_df['uniprot'])]


# print results of mapping
w_sifts = mapped_acc_df[(mapped_acc_df['blast']) & (mapped_acc_df['sifts'])].shape[0]
w_blast = mapped_acc_df[mapped_acc_df['blast']].shape[0]
tot = mapped_acc_df.shape[0]

print(" {} / {} ({:.2f} %) with blast results.".format(w_blast, tot, 100*w_blast/tot))
print(" {} / {} ({:.2f} %) with sifts results.".format(w_sifts, tot, 100*w_sifts/tot))



##################################
#### Import / Preprocess PDB  ####
##################################

# Import PDB Archive directory structure
pdbstore = clumpsptm.PdbStore(PDB_DIR) # NOTE: need to be absolutely sure there's no stray files in the directory-structure
print(pdbstore)

# Check for Missing PDBs and attempt to download
def check_missing_pdbs(sifts_filt_df):
    """Check how many PDBs are missing."""
    total_pdbs = len(set(sifts_filt_df['PDB']))
    missing_pdbs = set(sifts_filt_df['PDB']) - pdbstore.downloaded_pdbs
    print("Missing {} / {} ({:.2f} %) PDB files.".format(len(missing_pdbs), total_pdbs, 100*len(missing_pdbs)/total_pdbs))
    return missing_pdbs

missing_pdbs = check_missing_pdbs(sifts_filt_df)
# # attempt to download missing pdbs
# pdbstore.download_missing_pdbs(missing_pdbs)

# Drop PDBs unable to be downloaded
sifts_filt_df = sifts_filt_df[sifts_filt_df["PDB"].isin(pdbstore.downloaded_pdbs)]
if sifts_filt_df.shape[0]==0:
    raise Exception("No SIFTS hits had matching PDB files")

print("#### PDB IMPORT COMPLETE  --- ", datetime.now())



##################################
#### Parse Variable-Site File ####
##################################

def grab_ptm_site(x):
    """Grab actual PTM site."""
    res = list()
    for s in x.split(" "):
        if s.startswith(("K","Y","S","T")):
            res.append(s)
        else:
            pass
    if len(res) == 1:
        return res[0]
    else:
        return res

ptm_df = pmap_df.copy() # grab copy of pmap_df

# Expand PTMs in Each Feature
ptm_df[ptm_col] = [str.split(ptm_split) for str in ptm_df[ptm_col]] # convert variableSites column into list
#ptm_df[ptm_col] = ptm_df[ptm_col].apply(ast.literal_eval) # convert "site site site" format to {site, site, site} format by evaluating string as expression
ptm_df_long = ptm_df.explode(ptm_col).drop_duplicates() # make every list element into its own row


# Proteins with single PTM sites
ptm_sing_df = ptm_df_long.loc[ptm_df_long.index.map(lambda x: x.split("_")[-3]=="1" and x.split("_")[-4]=="1")].copy()
ptm_sing_df.loc[:,"ptmSite"] = ptm_sing_df[ptm_col].apply(grab_ptm_site)
# filter out malformed values
ptm_sing_df_filt = ptm_sing_df[ptm_sing_df['ptmSite'].notna()] # drop ptmSites that are NA
ptm_sing_df_filt = ptm_sing_df_filt[ptm_sing_df_filt['ptmSite'].apply(lambda x: len(x) > 0)] # drop ptmSites that are length zero (i.e. non K S T Y ptms)
if ptm_sing_df_filt.shape[0]==0:
	raise Exception("No single-site ptms in dataset")

if keep_multi_sites:
    # Proteins with multiple PTM sites
    ptm_multi_df = pd.concat((
        ptm_df_long.loc[ptm_df_long.index.map(lambda x: x.split("_")[-3]=="3" and x.split("_")[-4]=="3")].copy(),
        ptm_df_long.loc[ptm_df_long.index.map(lambda x: x.split("_")[-3]=="2" and x.split("_")[-4]=="2")].copy()
    ))
    ptm_multi_df.loc[:,"ptmSite"] = ptm_multi_df[ptm_col].apply(grab_ptm_site)
    ptm_multi_df = ptm_multi_df.explode("ptmSite")
    # filter out malformed values
    ptm_multi_df_filt = ptm_multi_df[ptm_multi_df['ptmSite'].notna()] # drop ptmSites that are NA
    ptm_multi_df_filt = ptm_multi_df_filt[ptm_multi_df_filt['ptmSite'].apply(lambda x: len(x) > 0)] # drop ptmSites that are length zero (i.e. non K S T Y ptms)
    if ptm_multi_df_filt.shape[0]==0:
    	warnings.warn("No multi-site ptms in dataset")
    # optionally filter sites that already exist in single-site data
    if filter_duplicate_sites:
        sing_keys = set(zip(ptm_sing_df_filt[accn_col], ptm_sing_df_filt["ptmSite"])) # get all unique PTM Sites in the single-site dataset
        ptm_multi_df_filt = ptm_multi_df_filt[ # filter out sites from the multisite dataset
            ~ptm_multi_df_filt[[accn_col,"ptmSite"]].apply(tuple, axis=1).isin(sing_keys) # that had a relevant single-site
        ]
    # Combine single-sites with multi-sites
    ptm_comb_df = pd.concat((ptm_sing_df_filt, ptm_multi_df_filt))
else:
    ptm_comb_df = ptm_sing_df_filt


ptm_comb_df = ptm_comb_df[ptm_comb_df[accn_col].isin(mapped_acc_df.index)] # subset to valid accession numbers
if ptm_comb_df.shape[0]==0:
	raise Exception("No overlap between database accession-numbers and BLASTed accession numbers")


# Join PTM features from source data-set to mapped accession-numbers
ptm_comb_df = pd.merge(
    ptm_comb_df.reset_index(), 
    mapped_acc_df.reset_index().rename(columns={"query":accn_col}),
    how="left",
    on=accn_col
).set_index(unique_id_col)


# ptm_comb_df["acc_res"] = ptm_comb_df["ptmSite"].apply(lambda x: x[0])
# ptm_comb_df["acc_res_i"] = ptm_comb_df["ptmSite"].apply(lambda x: int(x[1:-1]))
ptm_comb_df["acc_res"] = ptm_comb_df["ptmSite"].str.get(0)
ptm_comb_df["acc_res_i"] = ptm_comb_df["ptmSite"].str.slice(1,-1).astype(int)


# print out stats for your database
print("  * {} single PTM sites total in dataset".format(ptm_sing_df_filt.shape[0]))
if keep_multi_sites:
    print("  * {} multi PTM sites total in dataset".format(ptm_multi_df_filt.shape[0]))
else:
    print("  * Multi PTM sites will be excluded")


print("#### PTM PREPROCESSING COMPLETE  --- ", datetime.now())

##################################
#### Map Site IDX to Uniprot  ####
##################################

# Filter
ptm_comb_filt_df = ptm_comb_df[ptm_comb_df['Hsp_query-from'].notna()].copy()
ptm_comb_filt_df['Hsp_query-from'] = ptm_comb_filt_df['Hsp_query-from'].astype(int)
ptm_comb_filt_df['Hsp_hit-from'] = ptm_comb_filt_df['Hsp_hit-from'].astype(int)


def get_source_blast_idx(row):
    """
    Get source blast index.
    """
    d = dict()
    c = row['Hsp_query-from']
    for idx,res in enumerate(row['Hsp_qseq']):
        if res!='-':
            d[c] = idx
            c = c+1
    if row['acc_res_i'] not in d.keys():
        return "X"
    else:
        return d[row['acc_res_i']]

def get_uniprot_blast_i(row):
    """
    Get Uniprot Mapped Residue ID.
    """
    d = dict()
    c = row['Hsp_hit-from']
    for idx,res in enumerate(row['Hsp_hseq']):
        if res!='-':
            d[idx] = c
            c = c+1
    if row['acc_res_idx'] not in d.keys():
        return "X"
    else:
        return d[row['acc_res_idx']]

def _get_qres(row):
    """Get query residue"""
    try:
        return row["Hsp_qseq"][row["acc_res_idx"]]
    except:
        return None

def _get_hres(row):
    """Get query residue"""
    try:
        return row["Hsp_hseq"][row["acc_res_idx"]]
    except:
        return None

ptm_comb_filt_df["acc_res_idx"] = ptm_comb_filt_df.apply(get_source_blast_idx, 1)
ptm_comb_filt_df["acc_res"] = ptm_comb_filt_df.apply(_get_qres, 1)

# Get Uniprot Matches
ptm_comb_filt_df["uniprot_res"] = ptm_comb_filt_df.apply(_get_hres, 1)
ptm_comb_filt_df["uniprot_res_i"] = ptm_comb_filt_df.apply(get_uniprot_blast_i, 1)
ptm_comb_filt_df['uniprot_match'] = ptm_comb_filt_df['acc_res']==ptm_comb_filt_df['uniprot_res']



print("  * {} matched accession residues and uniprot residues".format(sum(ptm_comb_filt_df['uniprot_match'])))


if sum(ptm_comb_filt_df['uniprot_match'])==0:
    # raise Exception("No {} residues matched to Uniprot residues".format(accn_type))
    raise Exception("No residues matched to Uniprot residues")


print("#### PTM MAP-TO-UNIPROT COMPLETE  --- ", datetime.now())


################################################
#### Align Source Residues to PDB Residues  ####
################################################

sifts_filt_df['PDB_BEG'] = sifts_filt_df['PDB_BEG'].apply(lambda x: 0 if x=="None" else x) # overwrite "None" with 0
sifts_filt_df['len'] = sifts_filt_df['RES_END'] - sifts_filt_df['RES_BEG']
sifts_drop_dup = sifts_filt_df.drop_duplicates(subset=['PDB','len'])

ptm_comb_filt2_df = ptm_comb_filt_df[ptm_comb_filt_df["uniprot_match"]].copy()

def get_pdb_headers(sifts_df, pdbstore):
    """Get PDB headers.""" 
    pdb_chain_unis = [x for x in zip(sifts_df['PDB'], sifts_df['CHAIN'], sifts_df['SP_PRIMARY'])]
    print("   * Grabbing header info for {} pdb-chains".format(len(pdb_chain_unis)))
    headers = list()
    def get_header_info(pdb_ch_uni):
        """Download pdb."""
        pdb_headers = dict()
        try:
            pdb, chain, uniprot = pdb_ch_uni
            hd = pdbstore.load_header(pdb)
            for polymer in hd['polymers']:
                pdb_headers[polymer] = dict()
                if polymer.chid == chain:
                    for dbref in polymer.dbrefs:
                        pdb_headers[polymer]["PDB"] = pdb
                        pdb_headers[polymer]["CHAIN"] = chain
                        pdb_headers[polymer]["db_ref"] = dbref.database
                        pdb_headers[polymer]["db_accession"] = dbref.accession
                        pdb_headers[polymer]["db_first_from"] = dbref.first[0]
                        pdb_headers[polymer]["db_first_to"] = dbref.first[2]
                        if dbref.accession == uniprot:
                            continue
            headers.append(pd.DataFrame.from_dict(pdb_headers).dropna(axis=1))
        except:
            print("  * Error for {}".format(pdb_ch_uni))
    for pdb_chain_uni in tqdm(pdb_chain_unis):
        get_header_info(pdb_chain_uni)
    return headers

pdb_headers = pd.concat(get_pdb_headers(sifts_drop_dup, pdbstore),axis=1).T
pdb_headers.to_csv(os.path.join(REF_DIR, "pdb_headers.txt"), sep='\t')

pdb_headerpdb_headers = pd.read_csv(os.path.join(REF_DIR, "pdb_headers.txt"), sep='\t', index_col=0)

sifts_drop_dup = pd.merge(sifts_drop_dup.reset_index(), pdb_headers.reset_index().rename(columns={'index':'pdb_desc'}))
sifts_drop_dup['db_match'] = sifts_drop_dup['SP_PRIMARY']==sifts_drop_dup['db_accession']
sifts_drop_dup.groupby('db_match').size()



accession_numbers_to_use = np.intersect1d(
    mapped_acc_df[
        (mapped_acc_df['blast']) & (mapped_acc_df['sifts'])
    ].index,
    ptm_comb_filt2_df[accn_col]
)

if len(accession_numbers_to_use)==0:
	raise Exception("No valid accession numbers with BLAST hits, SIFTS hits, and matching residues.")


print("\n#### PDB HEADER COLLECTION COMPLETE  --- ", datetime.now())


##########################
#### Get PDB Matches  ####
##########################

# Code block where warnings are suppressed
# with warnings.catch_warnings():
#     warnings.simplefilter("ignore")
clumpsptm.mp.get_pdb_matches(
    accession_numbers_to_use,
    ptm_comb_filt2_df,
    sifts_drop_dup,
    pdbstore,
    os.path.join(REF_DIR, "pdb_matches"),
    n_threads,
    protein_id = accn_col
)
# NOTE: warnings are not an issue in version 3.7.3


print("#### PDB MATCHING COMPLETE  --- ", datetime.now())




############################
#### Filter & Finalize  ####
############################


# read in results
pdb_mapped_df = pd.concat([
    pd.read_parquet(x) for x in 
    glob.glob(os.path.join(REF_DIR, "pdb_matches/acc_mapped_sifts/*"))
])

# sort by percent_match
pdb_mapped_df['db_match'] = pdb_mapped_df['db_match'].astype(bool)
pdb_mapped_df = pdb_mapped_df[pdb_mapped_df['db_match']]
pdb_mapped_df = pdb_mapped_df[pdb_mapped_df['db_ref']=='UniProt']
pdb_mapped_df['percent_match'] = pdb_mapped_df['percent_match'].astype(float)
pdb_mapped_df = pdb_mapped_df.sort_values(['percent_match','len','proteins'], ascending=False)

# Select Top Structure by overlapping matched sites
pdb_mapped_filt_df = pdb_mapped_df.drop_duplicates(subset=['proteins'])

_prev = pdb_mapped_filt_df.shape[0]
pdb_mapped_filt_df = pdb_mapped_filt_df[pdb_mapped_filt_df['percent_match']>0]
pdb_mapped_filt_df['len'] = pdb_mapped_filt_df['len'].astype(int)
print(" {} / {} ({:.2f} %) structures with matching sites".format(
    pdb_mapped_filt_df.shape[0], _prev, 100*pdb_mapped_filt_df.shape[0]/_prev))




sites_mapped_df = list()

for idx,row in tqdm(pdb_mapped_filt_df.iterrows(), total=pdb_mapped_filt_df.shape[0]):
    _df = pd.read_parquet(
        os.path.join(REF_DIR, "pdb_matches/acc_mapped_sites", "{}.parquet".format(row["proteins"]))).loc[:,[
        "{}_{}_res_i".format(row["PDB"],row["CHAIN"]),
        "{}_{}_res".format(row["PDB"],row["CHAIN"]),
        "{}_{}_res_match".format(row["PDB"],row["CHAIN"])
    ]]
    _df.columns = ['pdb_res_i', 'pdb_res', 'pdb_res_match']
    _df['pdb'] = row["PDB"]
    _df['chain'] = row["CHAIN"]
    sites_mapped_df.append(_df)


sites_mapped_df = pd.concat(sites_mapped_df)

# Drop duplicates
sites_mapped_df = sites_mapped_df.reset_index().drop_duplicates().set_index(unique_id_col)

# Select for matched sites / rank multi sites by position
sites_mapped_df = sites_mapped_df[sites_mapped_df['pdb_res_match']]
sites_mapped_df['pdb_res_i_rank'] = sites_mapped_df.reset_index().groupby(unique_id_col)['pdb_res_i'].rank("dense").astype(int).values




# Drop duplicate ptmSites (i.e. NP_000005.2_K608k_1_1_608_608 - has two sites mapped from different cohorts)
ptm_comb_filt_pdb_df = ptm_comb_filt_df.reset_index().drop_duplicates([unique_id_col,'ptmSite']).set_index(unique_id_col)
# ptm_comb_filt_df.reset_index()[ptm_comb_filt_df.reset_index().duplicated(['id','ptmSite'], keep=False)] # see duplicated rows

# Rank multi-sites by position
ptm_comb_filt_pdb_df['pdb_res_i_rank'] = ptm_comb_filt_pdb_df.reset_index().groupby(unique_id_col)['acc_res_i'].rank("dense").astype(int).values

# Combine pdb level annotations
ptm_comb_filt_pdb_df = pd.merge(ptm_comb_filt_pdb_df.reset_index(), sites_mapped_df.reset_index(), how='left').set_index(unique_id_col)
ptm_comb_filt_pdb_df['pdb_res_match'] = ptm_comb_filt_pdb_df['pdb_res_match'].fillna(False)
ptm_comb_filt_pdb_df['pdb_res_match'] = (ptm_comb_filt_pdb_df['pdb_res_match']) & (ptm_comb_filt_pdb_df['pdb_res'] == ptm_comb_filt_pdb_df['acc_res'])

# Only matched sites
ptm_comb_filt_match_pdb_df = ptm_comb_filt_pdb_df[ptm_comb_filt_pdb_df['pdb_res_match']==True].copy()
ptm_comb_filt_match_pdb_df['pdb_res_i'] = ptm_comb_filt_match_pdb_df['pdb_res_i'].astype(int)



ptm_comb_filt_pdb_df.to_csv(os.path.join(OUT_DIR, output_prefix+"full_mapped_sites_to_pdbs.tsv"), sep='\t')
ptm_comb_filt_match_pdb_df.to_csv(os.path.join(OUT_DIR, output_prefix+"mapped_sites_to_pdbs.tsv"), sep='\t')


# add printout to end of file, so we know whether mapping was completed successfully
if (os.path.exists(os.path.join(OUT_DIR, output_prefix+"mapped_sites_to_pdbs.tsv")) and \
    os.path.exists(os.path.join(OUT_DIR, output_prefix+"full_mapped_sites_to_pdbs.tsv"))):
    print("INFO: Successfully wrote mapping files for Clumps PTM!")
else:
    raise Exception("Mapping files were not produced. Something has gone terribly wrong!")

# # copy file to local dir
# shutil.copyfile(os.path.join(REF_DIR, "full_mapped_sites_to_pdbs.tsv"), '/opt/input/full_mapped_sites_to_pdbs.tsv')
# shutil.copyfile(os.path.join(REF_DIR, "mapped_sites_to_pdbs.tsv"), '/opt/input/mapped_sites_to_pdbs.tsv')
