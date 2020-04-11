#!/bin/bash

for d in ../src/pgdac_*;
do
  mv "$d" "$(echo "$d" | sed s/pgdac/panoply/)";
done
