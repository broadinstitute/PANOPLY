#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

ns="broadcptacdev"

install_task_wdl(){
  task=$1
  hydrant install -m $task \
    -n $ns \
    -d tasks/$task/$task.wdl
}

install_wkflow_wdl(){
  wkflow=$1
  hydrant install -m $wkflow \
    -n $ns \
    -d workflows/$wkflow/$wkflow.wdl
}


# Uncomment for execution
# -------------------------------------

# install_task_wdl panoply_accumulate
# install_task_wdl panoply_association
# install_task_wdl panoply_cna_correlation
# install_task_wdl panoply_cna_setup
# install_task_wdl panoply_download
# install_task_wdl panoply_harmonize
# install_task_wdl panoply_normalize_ms_data
# install_task_wdl panoply_parse_sm_table
# install_task_wdl panoply_rna_protein_correlation
# install_task_wdl panoply_sampleqc
# install_task_wdl panoply_unified_pre


# install_wkflow_wdl panoply_main
