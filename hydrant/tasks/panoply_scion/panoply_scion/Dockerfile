FROM rocker/r-ver:4.2.1
MAINTAINER nclark@broadinstitute.org

# SCION requires R >= 4.0, but panoply_common is built on R 3.6.3, so this image is
# built from rocker/r-ver instead. Once panoply_common is upgraded, this can go back
# to being based on it (it already provides optparse, cmapR, etc.).

# System libraries
#   zlib1g-dev:       headers to compile cmapR's Bioconductor dependencies (Rhdf5lib, XVector)
#   libxml2, libglpk40: runtime libraries for igraph (a SCION dependency); libgmp10 comes with glpk
RUN apt-get update \
    && apt-get install -y --no-install-recommends zlib1g-dev libxml2 libglpk40 \
    && rm -rf /var/lib/apt/lists/*


# CRAN packages needed beyond SCION's own dependencies (installed with SCION below):
#   optparse (command-line parsing in the scion_*.R scripts), remotes (install_github),
#   BiocManager (Bioconductor repository URLs)
RUN install2.r --error --skipinstalled -n $(nproc) optparse remotes BiocManager \
    && rm -rf /tmp/downloaded_packages

# cmapR (reads the GCT inputs) is a Bioconductor package; its dependencies are compiled
# from source, so build in parallel. Bioconductor is pinned to 3.15 (June 2022): the default
# for R 4.2.1 is 3.16, whose cytolib needs a newer Boost (BH) than this CRAN snapshot has.
# BiocManager::install() won't pin an older release, so point install.packages() at the
# 3.15 repositories directly.
RUN Rscript -e 'options(Ncpus = parallel::detectCores()); \
      Sys.setenv(MAKEFLAGS = paste0("-j", parallel::detectCores())); \
      install.packages("cmapR", repos = BiocManager::repositories(version = "3.15"))' \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*


# SCION and its CRAN dependencies (cluster, doParallel, dtw, dtwclust, fastICA, ggplot2,
# igraph, randomForest, readr). Kept as the last R layer so changing scion_ref only
# rebuilds this step. MAKEFLAGS speeds up the one dependency (Matrix) that gets compiled
# from source here.
ARG scion_ref=master
RUN MAKEFLAGS="-j$(nproc)" Rscript -e "remotes::install_github('nmclark2/SCION', ref = '${scion_ref}', upgrade = 'never')" \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

COPY src /prot/proteomics/Projects/PGDAC/src
