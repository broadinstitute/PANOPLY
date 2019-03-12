#!/bin/bash

ver=5
name=`basename $PWD`
docker_tag=broadcptac/$name:$ver

mkdir src
#rm src/*
# copy code files from src-directory

cp -r ../../src/$name/* src
dos2unix src/*

# get ssGSEA scripts from GitHub
wget https://raw.githubusercontent.com/karstenkrug/ssGSEA2.0/master/src/ssGSEA2.0.R
dos2unix ssGSEA2.0.R
mv *.R src

wget https://raw.githubusercontent.com/karstenkrug/ssGSEA2.0/master/ssgsea-cli.R
dos2unix ssgsea-cli.R


docker build --rm -t $docker_tag .

#dos2unix *

rm -R src
rm *.R
