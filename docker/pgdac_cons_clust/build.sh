#!/bin/bash

ver=3
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

mkdir src

# copy code files from src-directory
cp ../../src/$name/* src
cp ../../src/assoc-analysis.r src
cp ../../src/run-pipeline.sh src

#cp -r ../../src/* src

dos2unix src/*

docker build --rm -t $docker_tag .

rm -R src/
