#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

echo -e "\n---------------------------"
echo -e "copy input files to Terra -"
echo -e "---------------------------\n"

source config.sh
source $src/tedmint-lib.sh

## check whether to split all or just one type
while getopts ":e:i:" opt; do
  case $opt in
    e) ext="$OPTARG";;
    i) input="$OPTARG";;
    \?) echo "Invalid Option -$OPTARG" >&2;;
  esac
done


datatype=`echo $input | rev | cut -d"." -f2- | cut -d"-" -f1 | rev`
bucket=$( get_bucket $wkspace $project )
gsutil -m cp split-data/$datatype/*.$ext gs://$bucket/$datatype/
