#!/bin/sh

# Commands to run to reasonably test the pipeline system
# run from test/results subdirectory

# 1. OPERATION inputSM
../../src/run-pipeline.sh inputSM -s ../input/proteome-medullo-v2-SMoutput-frag.ssv -e ../input/exptdesign.csv -r test-inputSM -c ../../src

# 2. OPERATION inputNorm
../../src/run-pipeline.sh inputNorm -n ../input/proteome-normalized-data.gct -r test-inputNorm -c ../../src

# 3. OPERATION RNAcorr 
../../src/run-pipeline.sh RNAcorr -f ../input/proteome-filtered-data.gct -r test-RNAcorr -c ../../src -rna ../input/microarray-data.gct
mv RNAcorr-output.tar RNAcorr-output1.tar
mv test-inputNorm test-inputNorm1
../../src/run-pipeline.sh RNAcorr -i inputNorm-output.tar -c ../../src -rna ../input/microarray-data.gct
