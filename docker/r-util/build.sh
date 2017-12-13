#!/bin/bash

ver=1
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

mkdir R-utilities
# cp from the mac results in fchmod permission denied errors -- ignore those and continue
cp -r /Volumes/prot_proteomics/LabMembers/manidr/R-utilities/*.[rR] R-utilities || true
cp -r /Volumes/prot_proteomics/LabMembers/manidr/R-utilities/GSEA R-utilities || true
docker build --rm -t $docker_tag .

rm -r R-utilities
