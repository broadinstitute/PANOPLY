#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#


function usage {
  # $1 SM table
  # $2 Expt design file
  # $3 Analysis directory
  # $4 Type
  # $5 RNA
  # $6 CNA
  echo "Usage: $0 <SM-table> <Expt-design-file> <Analysis-directory> <Type> <RNA-file> <CNA-file>"
}


if [[ $# -ne 6 ]]; then
  usage
  exit 1
fi


./common-code/run-pipeline.sh inputSM -s $1 -e $2 -r $3 -t $4 \
     -c common-code -d common-data -p config-custom.r -o tar/$3-sm.tar

./common-code/run-pipeline.sh normalize -i tar/$3-sm.tar  -t $4 \
     -c common-code -o tar/$3-norm.tar

./common-code/run-pipeline.sh RNAcorr -i tar/$3-norm.tar  -t $4 \
     -c common-code -rna $5 -o tar/$3-RNAcorr.tar

./common-code/run-pipeline.sh harmonize -i tar/$3-RNAcorr.tar  -t $4 \
     -c ./common-code -d ./common-data -rna $5 -cna $6 -o tar/$3-harm.tar

./common-code/run-pipeline.sh sampleQC -i tar/$3-harm.tar  -t $4 \
     -c common-code -o tar/$3-qc.tar

./common-code/run-pipeline.sh CNAsetup -i tar/$3-qc.tar  -t $4 \
     -c ./common-code -o tar/$3-CNAsetup.tar

./common-code/run-pipeline.sh CNAcorr -i tar/$3-CNAsetup.tar  -t $4 \
     -o tar/$3-CNAcorr.tar

./common-code/run-pipeline.sh assoc -i tar/$3-CNAcorr.tar  -t $4 \
     -c ./common-code -o tar/$3-assoc.tar

./common-code/run-pipeline.sh cluster -i tar/$3-assoc.tar  -t $4 \
     -c ./common-code -o tar/$3-final.tar
