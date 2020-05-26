FROM broadcptacdev/panoply_common:70df9e7
MAINTAINER rkothadi@broadinstitute.org

RUN Rscript -e "pacman::p_load ('prada')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('rhdf5')"
RUN Rscript -e "devtools::install_github('cmap/cmapR')"

COPY src /prot/proteomics/Projects/PGDAC/src
