#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
while getopts ":g:p:" opt; do
    case $opt in
        g) genome="$OPTARG";;
        p) proteome="$OPTARG";;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

mkdir genome && tar -xf $genome -C genome --strip-components 1
mkdir proteome && tar -xf $proteome -C proteome --strip-components 1
