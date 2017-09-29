#!/bin/bash

ver=3
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

mkdir src

cp -r ../../src/customProDB.r src
cp -r ../../src/aggregate_fasta.r src
cp -r ../../src/rmd-normalize.r src

docker build --rm -t $docker_tag .

rm -R src
