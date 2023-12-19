#!/bin/bash
set -e # exit upon error condition
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

echo -e "\n---------------------------"
echo -e "creates set memberships ---"
echo -e "---------------------------\n"

source config.sh

IFS='%' read -ra ssets <<< "$sets"
for sset in "${ssets[@]}"
do
  Rscript --verbose $src/r-source/create-sets.r \
    -s $sset \
    -l $dst \
    -w $wkspace
done
