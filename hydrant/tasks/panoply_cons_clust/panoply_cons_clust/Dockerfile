FROM broadcptacdev/panoply_common:74e21ac

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker for single-omics kmeans consensus clustering."

FROM broadcptacdev/panoply_common:latest
MAINTAINER rkothadi@broadinstitute.org

RUN rm -R /prot/proteomics/Projects/PGDAC/src/*
COPY src /prot/proteomics/Projects/PGDAC/src
