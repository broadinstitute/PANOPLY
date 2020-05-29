#!/bin/bash

echo -e "\n---------------------------"
echo -e "add links to input files --"
echo -e "---------------------------\n"

source config.sh
source $src/tedmint-lib.sh

bucket=$( get_bucket $wkspace $project )
Rscript --verbose $src/r-source/attribute-set.r \
  -b $bucket \
  -w $wkspace \
  -p $project \
  -c $csv_types \
  -g $gct_types \
  -o $dst
