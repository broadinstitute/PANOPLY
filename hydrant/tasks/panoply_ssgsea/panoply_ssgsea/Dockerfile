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

# additional CRAN packages needed:
#   optparse (command-line parsing), remotes (install_github), tidyverse
RUN install2.r --error --skipinstalled -n $(nproc) optparse remotes tidyverse

# cmapR (reads the GCT inputs) and ComplexHeatmap
# These are installed as prebuilt Linux binaries from Posit Package Manager
# with dependencies pulled both from frozen BioConductor and CRAN repositories
RUN Rscript -e 'options(Ncpus = parallel::detectCores()); \
      install.packages(c("cmapR", "ComplexHeatmap"), repos = c(getOption("repos"), "https://p3m.dev/bioconductor/__linux__/noble/2026-06-23/packages/3.23/bioc")); \
      stopifnot(requireNamespace("cmapR", quietly = TRUE), requireNamespace("ComplexHeatmap", quietly = TRUE))' # stop if packages didn't build

## install ssGSEA3.0 from GitHub (all of its dependencies are installed above or come from CRAN)
RUN Rscript -e 'remotes::install_github("broadinstitute/ssGSEA3.0", build_vignettes = FALSE); \
      stopifnot(requireNamespace("ssGSEA3", quietly = TRUE))' # stop if packages didn't build

## cleanup files from package installs
RUN rm -rf /tmp/downloaded_packages /tmp/Rtmp*