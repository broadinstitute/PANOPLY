# FROM broadcptacdev/panoply_common:latest
FROM rocker/r-ver:4.6.0
MAINTAINER proteogenomics@broadinstitute.org

LABEL description="Docker running gene-centric single-sample Gene Set Enrichment Analysis (ssGSEA) and PTM Signature Enrichment Analysis (PTM-SEA)".

ENV wd=/home/pgdac/

COPY . $wd
WORKDIR $wd

# ssGSEA3.0 requires R >= 4.1, but panoply_common is built on R 3.6.3, so this image is
# built from rocker/r-ver instead. Whenever you have time to fully revamp the panoply_common
# docker, this image can also be rebased...

# System libraries
#   zlib1g-dev:		headers to compile cmapR's Bioconductor dependencies (Rhdf5lib, XVector)
RUN apt-get update \
	&& apt-get install -y --no-install-recommends zlib1g-dev \
	&& rm -rf /var/lib/apt/lists/*

# additional CRAN packages needed:
#   optparse (command-line parsing), remotes (install_github),
#   BiocManager (Bioconductor repository URLs), tidyverse
RUN install2.r --error --skipinstalled -n $(nproc) optparse remotes BiocManager tidyverse \
    && rm -rf /tmp/downloaded_packages

# cmapR (reads the GCT inputs) from Bioconductor, pinned to repository v3.23 for stability
RUN Rscript -e 'options(Ncpus = parallel::detectCores()); \
      Sys.setenv(MAKEFLAGS = paste0("-j", parallel::detectCores())); \
      install.packages(c("cmapR", "ComplexHeatmap"), repos = BiocManager::repositories(version = "3.23"))' \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*


## install ssGSEA scripts from GitHub
RUN Rscript -e 'remotes::install_github("broadinstitute/ssGSEA3.0", build_vignettes = FALSE)'  