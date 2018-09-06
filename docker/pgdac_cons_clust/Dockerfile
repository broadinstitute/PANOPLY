FROM broadcptac/pgdac_cons_clust:2

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker for single-omics kmeans consensus clustering."


ENV wd /home/pgdac/
WORKDIR $wd

RUN rm -R /prot/proteomics/Projects/PGDAC/src/*
COPY src/* /prot/proteomics/Projects/PGDAC/src/

