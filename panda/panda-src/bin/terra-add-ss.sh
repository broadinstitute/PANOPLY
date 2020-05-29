#!/bin/bash

echo -e "\n---------------------------"
echo -e "Adding sample-set models --"
echo -e "---------------------------\n"

source config.sh
source $src/tedmint-lib.sh

## Add set names to config.sh
set_names=
set_names_str=""
IFS='%' read -ra ssets <<< "$sets"
for sset in "${ssets[@]}"
do
  IFS=',' read -ra temp <<< "$sset"
  set_names_str=$set_names_str"${temp[0]}"
  set_names_str=$set_names_str";"
  set_names+=( "${temp[0]}" )
done
echo -e "set_names=\""$set_names_str"\"" >> config.sh

Rscript --verbose $src/r-source/populate_sample_set.r \
  -g $gct_types \
  -c $csv_types \
  -s $dst/sample-set-members

mkdir -p sample-set-members-store
mv sample-set-members/* sample-set-members-store/.
fissfc entity_import -w $wkspace -p $project -f sample_set_membership.tsv
fissfc entity_import -w $wkspace -p $project -f sample_set_entity.tsv
