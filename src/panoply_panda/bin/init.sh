#!/bin/bash

echo -e "\n---------------------------"
echo -e "initializing w-space ------"
echo -e "---------------------------\n"

src=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
dst=`pwd`;
inp=$dst/pipeline-input;

## load libraries
source $src/tedmint-lib.sh
if [[ ! -f $dst/config.yaml ]]; then
  echo -e "$err config.yaml file not found. Exiting."
  exit
fi
Rscript --verbose $src/r-source/read-yaml.r -i config.yaml
source config.sh

## clean old config.sh and add params
sed '/## Others/,$d' $dst/config.sh > temp.sh
mv temp.sh config.sh
cat >> config.sh <<- EOM

## Others
src="$src"
dst="$dst"
inp="$inp"
EOM


## setup dir structure
if [[ -d "$inp" ]]; then
  rm -rf $inp;
fi

mkdir -p $inp
mkdir -p sample-set-members
mkdir -p additional-params
mkdir -p split-data
mkdir -p other-attributes
