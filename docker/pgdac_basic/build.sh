#!/bin/bash

ver=1
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

cp -r ../../src .
cp -r ../../data .
docker build --rm -t $docker_tag .

rm -rf src
rm -rf data
