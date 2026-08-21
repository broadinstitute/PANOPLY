FROM broadcptacdev/panoply_common:latest
MAINTAINER nclark@broadinstitute.org

# panoply_common already provides cluster, randomForest, doParallel, optparse,
# tidyverse (incl. ggplot2), devtools, and cmapR. Add what SCION needs beyond that.
RUN R -e "install.packages(c('dtwclust', 'fastICA', 'igraph', 'ranger'), \
  repos = 'https://cloud.r-project.org')"

ARG scion_ref=package-restructure
RUN R -e "devtools::install_github('nmclark2/SCION', ref = '${scion_ref}')"

COPY src /prot/proteomics/Projects/PGDAC/src
