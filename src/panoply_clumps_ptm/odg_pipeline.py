# load FASTA files
# os.system('pip install requests')

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

import inspect # for debugging
import shutil # for file copy


####################
####   Files    ####
####################
# Reference Files
UNIPROT_SWISSPROT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
SIFTS_DB = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"

# Source Files
INPUT_DIR = "/opt/input/"
CPTAC_REFSEQ_FASTA = os.path.join(INPUT_DIR,"Ensembl.human.hg19.clean3nr.602contams_20230913.fasta")
CPTAC_FEATURE_FILE = os.path.join(INPUT_DIR,"ODG_v3_pSTY_varsite.csv")

# Downloaded Database
PDB_DIR = "pdbs" # directory with PDB database

# Output Directory
REF_DIR = "example_ensembl"
fasta_folder = "ensembl_fasta"
os.makedirs(os.path.join(REF_DIR, fasta_folder), exist_ok=True)

# Other
n_threads = 12
ptm_split = " "
accn_col = "id.description"
accn_type = 'ENSEMBL'


###########################
####    File Import    ####
###########################
# Import Reference Files
clumpsptm.mp.dl_ref(REF_DIR, [UNIPROT_SWISSPROT, SIFTS_DB])

# Import PTM-Site Annotations & get list of accession numbers
pmap_df = pd.read_csv(CPTAC_FEATURE_FILE, index_col=0)
#pmap_df['accn_base'] = [re.sub("\\.\d+?","",accn) for accn in pmap_df['id.description']] # add base accession number to pmap_df
accn_arr = pmap_df[accn_col].drop_duplicates() # get unique accession numbers

# Import Sifts Data-Base
sifts_df = pd.read_csv(os.path.join(REF_DIR, "pdb_chain_uniprot.tsv"), comment="#", sep='\t', low_memory=False)


##################################
####   Read in FASTA Files    ####
##################################

# split FASTAs into individual files
clumpsptm.mp.split_fastas(
    CPTAC_REFSEQ_FASTA,
    os.path.join(REF_DIR, fasta_folder),
    naming="gencode" # determines which filename-splitting delimiter is used; "cptac"=" ", "gencode"="|"
)
individual_fastas_all = glob.glob(os.path.join(REF_DIR, fasta_folder, "*"))

# filter FASTAs to accession numbers in our dataset
#filt_in_df = [any(map(fasta.__contains__, accn_arr)) for fasta in individual_fastas_all]
pattern = re.compile('|'.join(map(re.escape, accn_arr)))
filt_in_df = [bool(pattern.search(fasta)) for fasta in individual_fastas_all] # slightly faster than map approach
#itertools.compress(individual_fastas_all,filt_in_df)
individual_fastas = list(itertools.compress(individual_fastas_all,filt_in_df))


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
# 		fn = "{dir}/{subdir}/{accn}.seq".format(accn=accn, dir=REF_DIR, subdir=fasta_folder)
# 		with open(fn, "w") as f:
# 			f.write(r.text)

# individual_fastas = glob.glob(os.path.join(REF_DIR, fasta_folder, "*"))
# #individual_fastas[:10]
# if len(individual_fastas) == 0:
# 	raise Exception("No FASTA files available for BLAST")
# os.listdir(os.path.join(REF_DIR, fasta_folder))

##################################
####    BLAST FASTA Files     ####
##################################
clumpsptm.mp.blast_sequences(
    os.path.join(REF_DIR, "uniprot_sprot.fasta"),
    individual_fastas,
    output_dir=REF_DIR,
    n_threads=n_threads,
    db_title="UniprotDB",
    blast_dir_name="refseq_to_uniprot_blast",
    seq_db_name="uniprot_db"
)

# collect BLASTed files
blasted_files_all = glob.glob(os.path.join(REF_DIR, "refseq_to_uniprot_blast/*"))
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
pdbstore = clumpsptm.PdbStore(PDB_DIR)
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
ptm_df["variableSites"] = [str.split(ptm_split) for str in ptm_df['variableSites']] # convert variableSites column into list
#ptm_df['variableSites'] = ptm_df['variableSites'].apply(ast.literal_eval) # convert "site site site" format to {site, site, site} format by evaluating string as expression
ptm_df_long = ptm_df.explode('variableSites').drop_duplicates() # make every list element into its own row


# Proteins with single PTM sites
ptm_sing_df = ptm_df_long.loc[ptm_df_long.index.map(lambda x: x.split("_")[-3]=="1" and x.split("_")[-4]=="1")].copy()
ptm_sing_df.loc[:,"ptmSite"] = ptm_sing_df['variableSites'].apply(grab_ptm_site)
# filter out malformed values
ptm_sing_df_filt = ptm_sing_df[ptm_sing_df['ptmSite'].notna()] # drop ptmSites that are NA
ptm_sing_df_filt = ptm_sing_df_filt[ptm_sing_df_filt['ptmSite'].apply(lambda x: len(x) > 0)] # drop ptmSites that are length zero (i.e. non K S T Y ptms)
if ptm_sing_df_filt.shape[0]==0:
	raise Exception("No single-site ptms in dataset")

# Proteins with multiple PTM sites
ptm_multi_df = pd.concat((
    ptm_df_long.loc[ptm_df_long.index.map(lambda x: x.split("_")[-3]=="3" and x.split("_")[-4]=="3")].copy(),
    ptm_df_long.loc[ptm_df_long.index.map(lambda x: x.split("_")[-3]=="2" and x.split("_")[-4]=="2")].copy()
))

ptm_multi_df.loc[:,"ptmSite"] = ptm_multi_df['variableSites'].apply(grab_ptm_site)
ptm_multi_df = ptm_multi_df.explode("ptmSite")

# filter out malformed values
ptm_multi_df_filt = ptm_multi_df[ptm_multi_df['ptmSite'].notna()] # drop ptmSites that are NA
ptm_multi_df_filt = ptm_multi_df_filt[ptm_multi_df_filt['ptmSite'].apply(lambda x: len(x) > 0)] # drop ptmSites that are length zero (i.e. non K S T Y ptms)
if ptm_multi_df_filt.shape[0]==0:
	raise Exception("No multi-site ptms in dataset")

# Combine & intersect accession numbers
ptm_comb_df = pd.concat((ptm_sing_df_filt, ptm_multi_df_filt))
ptm_comb_df = ptm_comb_df[ptm_comb_df[accn_col].isin(mapped_acc_df.index)] # subset to valid accession numbers
if ptm_comb_df.shape[0]==0:
	raise Exception("No overlap between database accession-numbers and BLASTed accession numbers")


# Join PTM features from source data-set to mapped accession-numbers
ptm_comb_df = pd.merge(
    ptm_comb_df.reset_index(), 
    mapped_acc_df.reset_index().rename(columns={"query":accn_col}),
    how="left",
    on=accn_col
).set_index("id")


# ptm_comb_df["acc_res"] = ptm_comb_df["ptmSite"].apply(lambda x: x[0])
# ptm_comb_df["acc_res_i"] = ptm_comb_df["ptmSite"].apply(lambda x: int(x[1:-1]))
ptm_comb_df["acc_res"] = ptm_comb_df["ptmSite"].str.get(0)
ptm_comb_df["acc_res_i"] = ptm_comb_df["ptmSite"].str.slice(1,-1).astype(int)


# print out stats for your database
print("  * {} single PTM sites total in dataset".format(ptm_sing_df_filt.shape[0]))
print("  * {} multi PTM sites total in dataset".format(ptm_multi_df_filt.shape[0]))


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
        if res is not '-':
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
        if res is not '-':
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

#print("  * {} matched accession residues and uniprot residues".format(sum(ptm_comb_filt_df['uniprot_match'])))


if sum(ptm_comb_filt_df['uniprot_match'])==0:
	raise Exception("No {} residues matched to Uniprot residues".format(accn_type))




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


##########################
#### Get PDB Matches  ####
##########################

clumpsptm.mp.get_pdb_matches(
    accession_numbers_to_use,
    ptm_comb_filt2_df,
    sifts_drop_dup,
    pdbstore,
    os.path.join(REF_DIR, "pdb_matches"),
    n_threads,
    protein_id = accn_col
)
# toDo: warnings in this are SUPER noisy. see if there's an easy fix you could push...






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
sites_mapped_df = sites_mapped_df.reset_index().drop_duplicates().set_index("id")

# Select for matched sites / rank multi sites by position
sites_mapped_df = sites_mapped_df[sites_mapped_df['pdb_res_match']]
sites_mapped_df['pdb_res_i_rank'] = sites_mapped_df.reset_index().groupby("id")['pdb_res_i'].rank("dense").astype(int).values




# Drop duplicate ptmSites (i.e. NP_000005.2_K608k_1_1_608_608 - has two sites mapped from different cohorts)
ptm_comb_filt_pdb_df = ptm_comb_filt_df.reset_index().drop_duplicates(['id','ptmSite']).set_index('id')
# ptm_comb_filt_df.reset_index()[ptm_comb_filt_df.reset_index().duplicated(['id','ptmSite'], keep=False)] # see duplicated rows

# Rank multi-sites by position
ptm_comb_filt_pdb_df['pdb_res_i_rank'] = ptm_comb_filt_pdb_df.reset_index().groupby("id")['acc_res_i'].rank("dense").astype(int).values

# Combine pdb level annotations
ptm_comb_filt_pdb_df = pd.merge(ptm_comb_filt_pdb_df.reset_index(), sites_mapped_df.reset_index(), how='left').set_index('id')
ptm_comb_filt_pdb_df['pdb_res_match'] = ptm_comb_filt_pdb_df['pdb_res_match'].fillna(False)
ptm_comb_filt_pdb_df['pdb_res_match'] = (ptm_comb_filt_pdb_df['pdb_res_match']) & (ptm_comb_filt_pdb_df['pdb_res'] == ptm_comb_filt_pdb_df['acc_res'])

# Only matched sites
ptm_comb_filt_match_pdb_df = ptm_comb_filt_pdb_df[ptm_comb_filt_pdb_df['pdb_res_match']==True].copy()
ptm_comb_filt_match_pdb_df['pdb_res_i'] = ptm_comb_filt_match_pdb_df['pdb_res_i'].astype(int)



ptm_comb_filt_pdb_df.to_csv(os.path.join(REF_DIR, "full_mapped_sites_to_pdbs.tsv"), sep='\t')
ptm_comb_filt_match_pdb_df.to_csv(os.path.join(REF_DIR, "mapped_sites_to_pdbs.tsv"), sep='\t')

# # copy file to local dir
# shutil.copyfile(os.path.join(REF_DIR, "full_mapped_sites_to_pdbs.tsv"), '/opt/input/full_mapped_sites_to_pdbs.tsv')
# shutil.copyfile(os.path.join(REF_DIR, "mapped_sites_to_pdbs.tsv"), '/opt/input/mapped_sites_to_pdbs.tsv')
