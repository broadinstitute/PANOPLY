FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN apt-get update --allow-releaseinfo-change
RUN echo Y | apt-get install libxml2

RUN R -e 'pacman::p_load ( Biostrings, seqinr, maftools, ActiveDriver )'

COPY src /prot/proteomics/Projects/PGDAC/src
