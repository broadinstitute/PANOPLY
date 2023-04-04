FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN apt-get update --allow-releaseinfo-change
RUN echo Y | apt-get install libxml2

RUN R -e 'pacman::p_load ( Biostrings, sequinr )'

RUN echo 'if (! require (mimp)) { \
        devtools::install_github("omarwagih/rmimp", INSTALL_opts="--no-staged-install")};' \
      > /tmp/mimp-packages.R \
   && Rscript /tmp/mimp-packages.R

COPY src /prot/proteomics/Projects/PGDAC/src
