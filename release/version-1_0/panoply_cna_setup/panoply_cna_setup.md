Documentation at https://github.com/broadinstitute/PANOPLY/blob/release-1_0/release/version-1_0/panoply_cna_setup/panoply_cna_setup.md

# ```panoply_cna_setup```

## Description

This module sets up directories and code for running CNA analysis.

## Input

Required inputs:

* ```tarball```: (`.tar` file) tarball from ```panoply_harmonize```
* ```type```: (String) data type
* ```yaml```: (`.yaml` file) master-parameters.yaml

Optional inputs:

* ```groupsFile```: (`.csv` file) subset of sample annotations
* ```peMaxDefault```: (Int, default = 250) default maximum processors/jobs for parallelism
* ```minCnaN```: (Int, default = 5) minimum number of samples required to perform CNA correlation analysis
* ```outFile```: (String, default = "panoply_cna_setup-output.tar") output `.tar` file name

## Output

Tarball including the following files in the `cna` subdirectory:

* Harmonized data matrices for CNA (`all-cna-matrix.csv`), RNA (`all-rna-matrix.csv`) and proteome (`all-pome-matrix.csv`) 
* A text file `file_table.tsv` listing all the data matrix files
