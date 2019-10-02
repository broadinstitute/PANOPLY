#!/bin/bash

git clone https://github.com/broadinstitute/pgdac_tedmint.git
mv pgdac_tedmint/src src
rm -rf pgdac_tedmint
docker build --rm -t broadcptac/pgdac_tedmint:1 .
rm -rf src
docker push broadcptac/pgdac_tedmint:1
