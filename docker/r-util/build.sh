#!/bin/bash

ver=1
name=`basename $PWD`
docker_tag=broadpgdac/$name:$ver

cp -r /prot/proteomics/Projects/R-utilities .
docker build --rm -t $docker_tag .

rm -r R-utilities
