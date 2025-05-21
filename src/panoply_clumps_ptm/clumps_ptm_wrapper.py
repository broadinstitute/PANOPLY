## Wrapper Method for clumpsptm function


##################################
####    Package Management    ####
##################################

# import Clumps PTM function
from clumpsptm.__main__ import main  # adjust path as needed

# parameters
import argparse
import yaml

# data wrangling
import pandas as pd

# utility
import pprint
import os
import sys

from datetime import datetime
# import random
import json

# configure ProDy to not printout messages, to avoid
from prody import confProDy
confProDy(verbosity='error')


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


### import command line parameters
parser = argparse.ArgumentParser(prog = 'ClumpsPTM Wrapper', description='Parameter-wrangling for automating CLUMPS-PTM.')
parser.add_argument("-y", "--yaml", type=str, help="Path to .yaml file with parameters.", required=True) # Parameter file with PANOPLY-defaults for Clumps-PTM

# file inputs
parser.add_argument('-i', "--diff_exp_file", '--input', type=str, required=True, help='<Required> Differential-expression input file (must have "id" column labelling comparisons).')
parser.add_argument('-m',"--var_sites_file", '--maps', type=str, required=True, help='<Required> Mapping between Site IDs and PDBs (index as indices that overlap input.')

# protein structures
parser.add_argument('-s','--pdbstore', type=str,required=True, help='<Required> path to PDBStore directory.')
# parser.add_argument('--alphafold', action='store_true', default=False, help='Run using alphafold structures.')
# parser.add_argument('--alphafold_threshold', type=float, default=75, help='Threshold confidence level for alphafold sites.')

# columnname parameters
parser.add_argument('-p', "--accession_col", '--protein_id', type=str, help='Unique protein id in input.')
parser.add_argument('-d', "--variable_sites_col", '--site_id', type=str, help='Unique site id in input.')
parser.add_argument('-w', "--weight_col", '--weight', type=str, help='<Required> Column with weights for CLUMPS-PTM in differential-expression file (ex. logFC).')

# other parameters
parser.add_argument('-c', '--run_combined', type=str2bool, help='Toggle for running all features combined.') 
parser.add_argument('-q', '--use_only_significant_sites', action='store_true', help='Only use significant sites for CLUMPS-PTM.')
parser.add_argument('--min_sites', type=str, help='Minimum number of sites.') # NOTE: placeholder-- currently cannot be toggled from commandline
# parser.add_argument('--min_sites', default=3, help='Minimum number of sites.')
parser.add_argument('-x', '--xpo', type=str, help='Soft distance threshold (t).') # NOTE: placeholder-- currently cannot be toggled from commandline
# parser.add_argument('-x', '--xpo', default=[6], type=list, help='Soft distance threshold (t).')
parser.add_argument('-n','--threads', type=str, default="1", help='Number of threads for sampling.') # NOTE: keep everything as string until it's passed to clumpsptm

# parser.add_argument('-r','--seed', type=str, default=None, help='Random seed to use in ClumpsPTM. WARNING: Only works for single-threaded processes.') 

# MANAGED BY WRAPPER
# parser.add_argument('-o','--output_dir', default=".", help='Output directory.')
# parser.add_argument('-f', '--features', nargs="*", default=None, help='Assays to subset for.') 
# parser.add_argument('-g', '--grouping', default=None, help='DE group to use.') 
# parser.add_argument('--subset', default=None, help='Subset sites.', choices=('positive','negative')) 

# testing / etc
parser.add_argument('-t','--test', action='store_true', default=False, help='Test run with n=5 proteins.')
parser.add_argument('-v', '--verbose', action='store_true', default=True, help='Verbosity.') # always run verbose (for now)

args = parser.parse_args()


# # testing arguments manually
# args = parser.parse_args([
#     "-y" "/opt/input/master-parameters.yaml",\
#     "-i" "/opt/input/filtered/full_de_cohort_cov_filt-to-1-4_filt-to-pdbs_filt-to-acKpSTY.tsv", \
#     "-m" "/opt/input/filtered/mapped_sites_to_pdbs_filt-to-pdbs.tsv", \
#     "--pdbstore", "pdbs", \
#     "--accession_col", "accession_number", \
#     "--weight_col", "gsea_rank", \
#     "-n", "12", \
#     # "-r", "2025", \
#     "-t"
# ])




### import default parameters from YAML
with open(args.yaml, 'r') as file:
    yaml_dict = yaml.safe_load(file)


# if (args.seed==None):  # SEED DOES NOT WORK
#     args.seed = datetime.now().strftime("%H%M%S") # set seed to current time
#     # args.seed = datetime.now().strftime("%Y%m%d%H%M%S") # seeds that are too long will fail


# override missing parameters with yaml defaults
if (args.accession_col==None):
    # args.accession_col = yaml_dict['global_parameters']['gene_mapping']['protein_id_col']
    args.accession_col = yaml_dict['panoply_ptm_normalization']['accession_number_colname']

if (args.variable_sites_col==None):
    args.variable_sites_col = yaml_dict['panoply_clumps_ptm']['mapping']['variable_sites_col']

if (args.weight_col==None):
    args.weight_col = yaml_dict['panoply_clumps_ptm']['analysis']['weight_col']

if (args.run_combined==None):
    args.run_combined = yaml_dict['panoply_clumps_ptm']['analysis']['run_combined']


# print parameters
print('\n\nParameters:')
pprint.pp(args.__dict__)
print('\n')

# write args dictionary to a JSON
os.makedirs('clumpsptm_runs')
with open(os.path.join("clumpsptm_runs",'params.json'), 'w') as f:
    json.dump(args.__dict__, f)


##################################
####      Data Wrangling      ####
##################################


## get all features and groups from diff_exp_file
_de = pd.read_csv(args.diff_exp_file, sep = '\t', low_memory=False)
groups = list(set(_de.id)) # get list of unique groupings
features = list(set(_de.feature)) # get list of unique features
if (args.run_combined):
    features.append('combined') # append 'combined' to run combined, if necessary


# Now call the actual CLI main function

for group in groups:
    for direction in ['positive', 'negative']:
        # Simulate CLI arguments (like from the shell)
        sys.argv = [
            'clumpsptm',  # dummy program name
            # file inputs
            '--input', str(args.diff_exp_file),
            '--maps', str(args.var_sites_file),
            # protein structures
            '--pdbstore', str(args.pdbstore),
            # columnname parameters
            '--protein_id', str(args.accession_col),
            '--site_id', str(args.variable_sites_col),
            '--weight',  str(args.weight_col),
            # groupings / subsets
            '--features', *features,
            '--grouping', str(group),
            '--subset', direction,
            # other parameters
            '--threads', str(args.threads),
            # '--seed', str(args.seed), # SEED ONLY WORKS FOR SINGLE-THREAD
            '--output_dir', os.path.join("clumpsptm_runs",str(group)+"_"+direction+"_results")
        ]
        # advanced
        if (args.use_only_significant_sites):
            sys.argv.append("-q")
        if (args.xpo!=None):
            sys.argv.append(["--xpo", args.xpo])
        if (args.min_sites!=None):
            sys.argv.append(["--min_sites", args.min_sites])
        # toggles
        if (args.test):
            sys.argv.append("-t")
        if (args.verbose):
            sys.argv.append("-v")
        # print command 
        print("## RUNNING COMMAND: \'"+' '.join(sys.argv)+"\'")
        # Run CLUMPS-PTM
        try:
            main()
            # subprocess.run(sys.argv)
        except ValueError as e: # if we get a value error
            if str(e) != 'NO RESULTS FILES FOUND.':
                raise
            else:
                print("## WARNING: NO RESULTS FILES FOUND FOR GROUP '"+str(group)+"'") # print a warning and move on


