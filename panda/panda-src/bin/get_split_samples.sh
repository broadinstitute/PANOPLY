#!/bin/bash

echo -e "\n---------------------------"
echo -e "get split samples ---------"
echo -e "---------------------------\n"

source config.sh
source $src/tedmint-lib.sh


fissfc sample_list -w $wkspace -p $project -t sample > samples_list.txt
declare -a samples
samples=(`cat "samples_list.txt"`)

> filled_samples_0.tsv
> filled_samples.tsv
for s in "${samples[@]}"
do 
  fissfc attr_get \
      -w $wkspace \
      -p $project \
      -t sample \
      -e $s >> filled_samples_0.tsv
done

head -n 1 filled_samples_0.tsv > filled_samples.tsv
sed '/^entity:sample_id/d' < filled_samples_0.tsv >> filled_samples.tsv
rm samples_list.txt filled_samples_0.tsv
