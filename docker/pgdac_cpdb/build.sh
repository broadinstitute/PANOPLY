#!/bin/bash

ver=3
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

mkdir src

cp -r ../../src/customProDB.r src
cp -r ../../src/aggregate_fasta.r src

docker build --rm -t $docker_tag .

rm -rf src
