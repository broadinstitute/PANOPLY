#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
while getopts ":i:o:r:m:" opt; do
    case $opt in
        i) tarball="$OPTARG";;
        o) output_tar="$OPTARG";;
        r) analysisDir="$OPTARG";;
        m) module="$OPTARG";;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

src=`pwd`
mkdir $analysisDir && tar -xvf $tarball -C $analysisDir --strip-components 1
cp /prot/proteomics/Projects/PGDAC/src/contrast-to-gct.r $analysisDir/${module}/.
cd $analysisDir/${module}
R CMD BATCH contrast-to-gct.r

tar -cvf $src/$output_tar contrasts/*
