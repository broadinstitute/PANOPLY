#!/bin/bash

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


## check files in bucket, and copy only those missing
datatype=`echo $input | rev | cut -d"." -f2- | cut -d"-" -f1 | rev`
bucket=$( get_bucket $wkspace $project )
bucket_files="tmp-$datatype-bucketlist.txt"
echo "" > $bucket_files   # create an empty file to store files missing in bucket
for f in `gsutil ls gs://$bucket/$datatype/*.$ext`; do basename $f >> $bucket_files; done 
missing=`comm -23 <( (basename -a split-data/$datatype/*.$ext | sort) ) <(sort $bucket_files)`
if [ "$missing" != "" ]; then
  for m in $missing; do
    gsutil -m cp split-data/$datatype/$m gs://$bucket/$datatype/$m
  done
fi
rm $bucket_files
