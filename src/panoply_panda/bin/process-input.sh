#!/bin/bash

echo -e "\n---------------------------"
echo -e "reading config file -------"
echo -e "---------------------------\n"

source config.sh

## create symbolic links
IFS=';' read -ra gct_files_a <<< "$gct_files"
IFS=';' read -ra gct_types_a <<< "$gct_types"
for idx in "${!gct_files_a[@]}"
do
  add="$inp/$wkspace-${gct_types_a[idx]}.gct"
  if [[ -z $gcts ]]; then
    gcts=$add
  else
    gcts="$gcts;$add"
  fi
  ln -s $freeze_path/${gct_files_a[idx]} $add
done

IFS=';' read -ra csv_files_a <<< "$csv_files"
IFS=';' read -ra csv_types_a <<< "$csv_types"
for idx in "${!csv_files_a[@]}"
do
  add="$inp/$wkspace-${csv_types_a[idx]}.csv"
  if [[ -z $csvs ]]; then
    csvs=$add
  else
    csvs="$csvs;$add"
  fi
  ln -s $freeze_path/${csv_files_a[idx]} $add
done


## add params to config
echo -e "gcts=\"$gcts\"" >> config.sh
echo -e "csvs=\"$csvs\"" >> config.sh
