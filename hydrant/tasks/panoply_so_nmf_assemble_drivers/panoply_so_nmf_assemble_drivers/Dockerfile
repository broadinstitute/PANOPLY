FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

LABEL description="Docker to concatenate nmf driver-genes into gct"

RUN apt-get update --allow-releaseinfo-change && apt-get install -y dos2unix

ENV wd /home/pgdac/
COPY . $wd
WORKDIR $wd
RUN dos2unix src/*
