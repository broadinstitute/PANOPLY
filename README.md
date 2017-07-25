# 
# PGDAC code and pipelines
#

The ```docker``` directory contains code to create required docker images for use with the PGDAC pipeline. 
The ```r-util``` image includes R code from https://github.com/broadinstitute/proteomics-Rutil, 
and is the basis for the ```broadcptac/pgdac_basic``` image.

The essential code for the ```pgdac_basic``` image is contained in the ```src``` directory, 
with the associated workflow in the ```wdl``` directory. ```firecloud``` contains 
FireCloud documentation and supporting files needed to run workflows.

In order to follow the instructions in ```docs/Firecloud workflows from R code modules.ipynb```, 
clone the ```gdac-firecloud``` repository from https://github.com/broadinstitute/gdac-firecloud
into the root ```PGDAC``` directory.

The ```pipeline-common-code``` is a precursor to ```src``` where code for the pipeline is developed
using actual projects as use-cases. This directory is maintained in order have a version of the 
pipeline code that can be run on the commandline, without FireCloud support.

Email manidr@broadinstitute.org with questions.
