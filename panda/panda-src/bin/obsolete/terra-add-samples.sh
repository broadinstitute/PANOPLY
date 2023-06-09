#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

echo -e "\n---------------------------"
echo -e "add samples to terra ------"
echo -e "---------------------------\n"

source config.sh
source $src/tedmint-lib.sh

Rscript --verbose $src/r-source/populate_sample_meta.r \
  -g $gct_types \
  -c $csv_types \
  -o $dst \
  -w $wkspace

bucket=$( get_bucket $wkspace $project )
fissfc entity_import -w $wkspace -p $project -f $dst/participant.tsv
fissfc entity_import -w $wkspace -p $project -f $dst/samples.tsv
