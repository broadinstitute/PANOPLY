#!/bin/bash

ver=1
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

cp -r ../../src .
docker build --rm -t $docker_tag .

rm -rf src
