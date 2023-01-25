#!/bin/bash

#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#


# get the file basenames, remove extension (only csv of gct allowed)
d1_file_name=$(basename -s .csv $(basename -s .gct $d1_file))
d2_file_name=$(basename -s .csv $(basename -s .gct $d2_file))
sample_file_name=$(basename -s .csv $(basename -s .gct $sample_file))

# write file with d1 and d2 file names to use as task outputs
echo $d1_file_name > d1_file_name.txt
echo $d2_file_name > d2_file_name.txt

# preprocessed data directory and names
preprocessed_dir="./cosmo_preprocessed_data"
d1_file_preprocessed="$preprocessed_dir/$d1_file_name.tsv"
d2_file_preprocessed="$preprocessed_dir/$d2_file_name.tsv"
sample_file_preprocessed="$preprocessed_dir/$sample_file_name.tsv"
sample_label_preprocessed="$preprocessed_dir/sample_label.txt"

# directory where the outputs will go
out_dir="./output"

# directory and files where the data to use in COSMO will go
data_use_dir="./data_use"
d1_file_use="$data_use_dir/$d1_file_name.tsv"
d2_file_use="$data_use_dir/$d2_file_name.tsv"
sample_file_use="$data_use_dir/$sample_file_name.tsv"

# folders for each method
method1_out_folder="$out_dir/method1_folder"
method2_out_folder="$out_dir/method2_folder"
final_res_out_folder="$out_dir/final_res_folder"

# make the appropriate folders
mkdir $out_dir
mkdir $data_use_dir

# run panoply-specific preprocessing
Rscript /prot/proteomics/Projects/PGDAC/src/cosmo/panoply_cosmo_preprocessing.R \
  $d1_file \
  $d2_file \
  $sample_file \
  $yaml_file

# load the sample labels
sample_label=$(cat $sample_label_preprocessed)

# cosmo format data
R -e \
	"source('/prot/proteomics/Projects/PGDAC/src/cosmo/cosmo/tools.R'); format_input_data('$d1_file_preprocessed', '$d2_file_preprocessed', '$sample_file_preprocessed', out_dir = '$data_use_dir')"

# METHOD 1
mkdir $method1_out_folder

R -e \
    "source('/prot/proteomics/Projects/PGDAC/src/cosmo/cosmo/method1_function.R');
    d1_file <- '$d1_file_use'; 
    d2_file <- '$d2_file_use'; 
    sample_file <- '$sample_file_use';
    gene_file <- '/prot/proteomics/Projects/PGDAC/src/cosmo/cosmo/genes.tsv';
    out_dir <- '$method1_out_folder';
    clinical_attributes <- unlist(strsplit(x='$sample_label',split=','));
    run_2b(d1_file, d2_file, sample_file, gene_file, out_dir=out_dir, clinical_attributes=clinical_attributes)"

# METHOD 2
mkdir $method2_out_folder

python /prot/proteomics/Projects/PGDAC/src/cosmo/cosmo/method2_function.py \
    -d1 $d1_file_use \
    -d2 $d2_file_use \
    -s $sample_file_use \
    -l $sample_label \
    -o $method2_out_folder

# COMBINE METHODS
mkdir $final_res_out_folder

R -e \
    "source('/prot/proteomics/Projects/PGDAC/src/cosmo/cosmo/method1_function.R');
    source('/prot/proteomics/Projects/PGDAC/src/cosmo/cosmo/combine_methods.R');
    method1_folder <- '$method1_out_folder';
    method2_folder <- '$method2_out_folder';
    out_dir <- '$final_res_out_folder';
    sample_annotation_file <- '$sample_file_use';
    clinical_attributes <- unlist(strsplit(x='$sample_label',split=','));
    combine_methods(method1_folder, method2_folder, sample_annotation_file, clinical_attributes = clinical_attributes, out_dir = out_dir, prefix = 'cosmo')"
