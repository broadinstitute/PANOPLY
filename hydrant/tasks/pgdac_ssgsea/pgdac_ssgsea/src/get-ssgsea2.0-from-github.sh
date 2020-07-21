#!/bin/bash

rm -R src
mkdir src

wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/ssGSEA2.0.R
dos2unix ssGSEA2.0.R
mv ssGSEA2.0.R src

wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/parse_yaml.R
dos2unix parse_yaml.R
mv parse_yaml.R src

wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/ssgsea-cli.R
dos2unix ssgsea-cli.R
