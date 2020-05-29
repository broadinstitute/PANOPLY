#!/bin/bash

cd
home=`pwd`
mkdir -p .tedmint
( 
  cd .tedmint;
  git clone https://github.com/broadinstitute/pgdac_tedmint.git 
)
