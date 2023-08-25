#!/bin/bash
set -eu
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

dir_name="ssgsea_assoc"
tar_array=("$@") # read all arguments as an array
for index in "${!tar_array[@]}" # for each index
do
  # make a temporary-index directory and untar file
  mkdir -p $dir_name/$index && tar -xf "${tar_array[$index]}" -C "$dir_name/$index";
  # get the group from the relevant parameters.txt file 
  group="$( echo -e "$( grep "input gct" $dir_name/$index/*parameters.txt | \
    cut -d ':' -f 2 | tr -d '[:space:]')" | \
    rev | cut -d '.' -f2- | rev )"
  # make a new directory with the group label, and remove tmp-index directory
  mkdir -p $dir_name/"ssgsea-$group";
  cp -r $dir_name/$index/* $dir_name/"ssgsea-$group"/.;
  rm -rf $dir_name/$index;
done

tar -cf ssgsea_assoc.tar ssgsea_assoc/