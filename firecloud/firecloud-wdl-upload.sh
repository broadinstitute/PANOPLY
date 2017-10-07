#!/bin/bash

# run firecloud docker to access api
docker run --rm -it -v "$HOME"/.config:/.config -v "$PWD":/working broadinstitute/firecloud-cli bash

# authenticate in 
gcloud auth login

# upload various workflows (wdl)
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_basic -t Workflow -y "PGDAC basic pipeline" pgdac_basic/pgdac_basic.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n cna_analysis -t Workflow -y "CNA-RNA-Proteome correlation analysis" cna_analysis/cna_analysis.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n cna_analysis_bysubgroup -t Workflow -y "CNA-RNA-Proteome correlaion by subgroup" cna_analysis_bysubgroup.wdl
