#!/bin/bash

ver=2
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

# cp from the mac results in fchmod permission denied errors -- ignore those and continue
cp -r ../../src . || true
rm -rf src/.Rproj* src/*.Rproj
dos2unix src/*

cp -r -L ../../data . || true
cp -r -L ../packages . || true
# in the above dereference any symbolic links using -L

# also update R-utilities to include changes without having to rebuild r-util
mkdir R-utilities
cp ../../../R-utilities/*.[rR] R-utilities || true
cp -r ../../../R-utilities/GSEA R-utilities || true

docker build --rm -t $docker_tag .

rm -rf src
rm -rf data
rm -rf packages
rm -rf R-utilities
