FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN Rscript -e "pacman::p_load (prada, cmapR)"

COPY src /prot/proteomics/Projects/PGDAC/src
