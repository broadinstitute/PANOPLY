#!/bin/bash

ver=1
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

mkdir R-utilities
cp -r /Volumes/prot_proteomics/LabMembers/manidr/R-utilities/*.[rR] R-utilities
cp -r /Volumes/prot_proteomics/LabMembers/manidr/R-utilities/GSEA R-utilities
docker build --rm -t $docker_tag .

rm -r R-utilities
