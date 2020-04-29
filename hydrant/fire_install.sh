#!/bin/bash
ns="broadcptac"
hydrant install -m panoply_main -n $ns -d workflows/panoply_main/panoply_main.wdl
#hydrant install -m panoply_download -n $ns -d tasks/panoply_download/panoply_download.wdl
#hydrant install -m panoply_parse_sm_table -n broadcptac -d tasks/panoply_parse_sm_table/panoply_parse_sm_table.wdl
#hydrant install -m panoply_normalize_ms_data -n $ns -d tasks/panoply_normalize_ms_data/panoply_normalize_ms_data.wdl
#hydrant install -m panoply_rna_protein_correlation -n $ns -d tasks/panoply_rna_protein_correlation/panoply_rna_protein_correlation.wdl
#hydrant install -m panoply_harmonize -n $ns -d tasks/panoply_harmonize/panoply_harmonize.wdl
#hydrant install -m panoply_sampleqc -n $ns -d tasks/panoply_sampleqc/panoply_sampleqc.wdl
#hydrant install -m panoply_cna_setup -n $ns -d tasks/panoply_cna_setup/panoply_cna_setup.wdl
#hydrant install -m panoply_cna_correlation -n $ns -d tasks/panoply_cna_correlation/panoply_cna_correlation.wdl
#hydrant install -m panoply_association -n $ns -d tasks/panoply_association/panoply_association.wdl

