FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

LABEL description="Docker to create tar ball as input for panoply_mo_nmf"

RUN apt-get update && apt-get install -y dos2unix

ENV wd /home/pgdac/
COPY . $wd
WORKDIR $wd
RUN dos2unix src/*
