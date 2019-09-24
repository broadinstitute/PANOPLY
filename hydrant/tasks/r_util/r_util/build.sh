#!/bin/bash

ver=2
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

rm -rf proteomics-Rutil
rm -rf R-utilities
git clone https://github.com/broadinstitute/proteomics-Rutil.git
mv proteomics-Rutil R-utilities

docker build --rm -t broadcptac/r_util:2 .

rm -r R-utilities
