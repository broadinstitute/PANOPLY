#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

ver=5
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

mkdir src

cp -r ../../src/customProDB.r src
cp -r ../../src/aggregate_fasta.r src

cp -r ../../src/rmd-normalize.r src
cp -r ../../src/rmd-rna-seq-correlation.r src

cp -r ../../src/ssGSEA2.0.R src
cp -r ../../src/run_ssGSEA2.0.R src

docker build --rm -t $docker_tag .

rm -R src
