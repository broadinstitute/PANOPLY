#!/bin/bash
ver=3

name=`basename $PWD`

docker_tag=broadcptac/$name:$ver

mkdir src

# copy code files from src-directory
cp -r ../../src/$name/* src

docker build --rm -t $docker_tag .

rm -R src