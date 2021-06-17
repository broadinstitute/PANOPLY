#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
while getopts ":a:s:" opt; do
    case $opt in
	      a) analysis_dir="$OPTARG";;
        s) ssgsea_assoc="$OPTARG";;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

scatter_processing()
{
  dir_name=$1
  array=($(ls $dir_name/*.tar))
  for index in "${!array[@]}"
  do
    mkdir -p $dir_name/$index && tar xf ${array[index]} -C $dir_name/$index;
    group="$( echo -e "$( grep "input gct" $dir_name/$index/*parameters.txt | \
      cut -d ':' -f 2 | tr -d '[:space:]')" | \
      rev | cut -d '.' -f2- | rev )"
    mkdir -p $dir_name/"ssgsea-$group";
    cp -r $dir_name/$index/* $dir_name/"ssgsea-$group"/.;
    rm -rf $dir_name/$index;
  done
}

dir_create()
{
  cd $src;
  scatter_processing $src/$ssgsea_assoc
}


collect()
{
  cd $src;
  
  mkdir -p $full_path;
  cp -r $ssgsea_assoc $full_path/$ssgsea_assoc/;

  rm -rf $ssgsea_assoc;

}

src=`pwd`

full="$analysis_dir-full-results"

full_path=$src/$full
dir_create
collect
cd $src;
tar -cvf panoply_main_full.tar $full