#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

echo -e "\n---------------------------"
echo -e "splitting input data-------"
echo -e "---------------------------\n"

## load libraries
source config.sh
source $src/tedmint-lib.sh

while getopts ":e:i:" opt; do
  case $opt in
    e) ext="$OPTARG";;
    i) input="$OPTARG";;
    \?) echo "Invalid Option -$OPTARG" >&2;;
  esac
done

datatype=`echo $input | rev | cut -d"." -f2- | cut -d"-" -f1 | rev`
mkdir -p split-data/$datatype
Rscript --verbose $src/r-source/split.r \
  -e $ext \
  -i $input \
  -d $datatype \
  -o split-data/$datatype
