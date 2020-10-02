#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

cd
home=`pwd`
mkdir -p .tedmint
( 
  cd .tedmint;
  git clone https://github.com/broadinstitute/pgdac_tedmint.git 
)
