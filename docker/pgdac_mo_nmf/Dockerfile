FROM rocker/r-ver:3.6.3

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker for PANOPLY's multi-omic NMF clustering."

ENV wd /home/pgdac/

COPY . $wd
WORKDIR $wd

RUN apt-get update
RUN apt-get install -y apt-utils
## required for R-pacakge 'XML'
RUN apt-get install -y libxml2-dev

## rmarkdown
RUN apt-get install -y libc6-dev
RUN apt-get install -y perl
RUN apt-get install -y pandoc

## ssl
RUN apt-get install -y libssl-dev
## libcurl
RUN apt-get install -y libcurl4-openssl-dev

RUN apt-get install zlib1g-dev

## install R-packages
RUN R -e "install.packages('pacman')"

## Bioconductor
RUN R -e "pacman::p_install('BiocManager')"

## devtools
RUN R -e "pacman::p_install('usethis')"
RUN R -e "pacman::p_install('git2r')"
RUN R -e "pacman::p_install('devtools')"

## cmapR
Run R -e "pacman::p_install('rhdf5')"
RUN R -e "devtools::install_github('cmap/cmapR')"

## stats
RUN R -e "pacman::p_install('limma')"
RUN R -e "pacman::p_install('statmod')"
RUN R -e "pacman::p_install('impute')"
RUN R -e "pacman::p_install('Rtsne')"

## parallel loops
RUN R -e "pacman::p_install('foreach')"
RUN R -e "pacman::p_install('doParallel')"
RUN R -e "pacman::p_install('magrittr')"

## tidyverse
RUN R -e "pacman::p_install('glue')"
RUN R -e "pacman::p_install('tidyselect')"
RUN R -e "pacman::p_install('dplyr')"
RUN R -e "pacman::p_install('tibble')"

## misc
RUN R -e "pacman::p_install('stringi')"
RUN R -e "pacman::p_install('pkgmaker')"
RUN R -e "pacman::p_install('reshape2')"
RUN R -e "pacman::p_install('yaml')"
RUN R -e "pacman::p_install('optparse')"
RUN R -e "pacman::p_install('verification')"
RUN R -e "pacman::p_install('gtools')"
RUN R -e "pacman::p_install('rngtools')"
RUN R -e "pacman::p_install('rmarkdown')"
RUN R -e "pacman::p_install('maptools')"

## ploting
RUN R -e "pacman::p_install('ggplot2')"
RUN R -e "pacman::p_install('plotly')"
RUN R -e "pacman::p_install('vioplot')"
RUN R -e "pacman::p_install('UpSetR')"
RUN R -e "pacman::p_install('ComplexHeatmap')"
RUN R -e "pacman::p_install('ggpubr')"
RUN R -e "install.packages('circlize')"
RUN R -e "pacman::p_install('pheatmap')"
RUN R -e "pacman::p_install('RColorBrewer')"

## NMF
RUN R -e "pacman::p_install('NMF')"
RUN R -e "pacman::p_install('ccfindR')"

## annotation packages
RUN R -e "pacman::p_install('IRanges')"
RUN R -e "pacman::p_install('S4Vectors')"
RUN R -e "pacman::p_install('annotationDbi')"
RUN R -e "pacman::p_install('org.Hs.eg.db')"
RUN R -e "pacman::p_install('org.Mm.eg.db')"
RUN R -e "pacman::p_install('org.Rn.eg.db')"

## I/O
RUN R -e "pacman::p_install('WriteXLS')"
RUN R -e "pacman::p_install('readxl')"
