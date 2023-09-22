#!/bin/bash
set -e # exit upon error condition
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

## load parameters
source config.sh

echo -e "\n---------------------------"
echo -e "creating groups file ------"
echo -e "---------------------------\n"

## create groups file from the specified columns
if [[ ! -z $groups_cols ]]; then
  Rscript --verbose $src/r-source/create-groups.r \
    -o $dst \
    -w $wkspace \
    -c $groups_cols
  echo -e "csv_types=\"$csv_types;groups\"" >> config.sh
  echo -e "csv_files=\"$csv_files;$wkspace-groups.csv\"" >> config.sh
  echo -e "csvs=\"$csvs;$inp/$wkspace-groups.csv\"" >> config.sh
fi
