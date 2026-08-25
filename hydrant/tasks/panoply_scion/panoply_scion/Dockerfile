FROM broadcptacdev/panoply_common:latest
MAINTAINER nclark@broadinstitute.org

# panoply_common already provides cluster, randomForest, doParallel, optparse,
# tidyverse (incl. ggplot2 and readr), devtools, and cmapR. Add what SCION
# needs beyond that. (ranger was dropped -- SCION's randomForest engine is
# the only one now; see nmclark2/SCION's v5.0 changelog.)
RUN R -e "install.packages(c('dtwclust', 'fastICA', 'igraph'), \
  repos = 'https://cloud.r-project.org')"

ARG scion_ref=master
RUN R -e "devtools::install_github('nmclark2/SCION', ref = '${scion_ref}')"

COPY src /prot/proteomics/Projects/PGDAC/src
