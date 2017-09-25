#!/bin/bash

ver=3
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

cp -r ../../src/customProDB.r .
cp -r ../../src/aggregate_fasta.r .

docker build --rm -t $docker_tag .

rm -rf src
