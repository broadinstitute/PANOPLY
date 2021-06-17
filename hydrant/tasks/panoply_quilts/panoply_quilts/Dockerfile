FROM python:3
MAINTAINER proteogenomics@broadinstitute.org

# get the pyQUILTS GitHub repo:
RUN git clone https://github.com/ekawaler/pyQUILTS.git /src/QUILTS/pyQUILTS/

# make the results directory
RUN mkdir /src/QUILTS/output/

# compile the read_chr_bed.c file:
RUN gcc -o /src/QUILTS/pyQUILTS/read_chr_bed /src/QUILTS/pyQUILTS/read_chr_bed.c

# cp the quilts_references.sh
COPY src /src/QUILTS/src
RUN chmod 777 /src/QUILTS/src/quilts_references.sh
