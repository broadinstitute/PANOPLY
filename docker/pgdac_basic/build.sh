#!/bin/bash

ver=1
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

# cp from the mac results in fchmod permission denied errors -- ignore those and continue
cp -r ../../src . || true
rm -rf src/.Rproj* src/*.Rproj
cp -r ../../data . || true
docker build --rm -t $docker_tag .

rm -rf src
rm -rf data
