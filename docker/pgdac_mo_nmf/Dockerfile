#FROM broadcptac/r-util:1
FROM r-base

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker for multi-omic NMF clustering."

ENV wd /home/pgdac/

COPY . $wd
#COPY src src/

WORKDIR $wd

RUN apt-get update
## required for R-pacakge 'XML'
RUN apt-get -t unstable install -y libxml2-dev
## rmarkdown
RUN apt-get install -y libc6-dev
RUN apt-get install -y perl
RUN apt-get -t unstable install -y pandoc
## ssl
RUN apt-get install -y libssl-dev
## libcurl
RUN apt-get install -y libcurl4-openssl-dev


## install R-packages
RUN R -e "install.packages('pacman')"
RUN R -e "library(pacman);p_install('stringi')"
RUN R -e "library(pacman);p_install('pkgmaker')"
RUN R -e "library(pacman);p_install('reshape2')"
RUN R -e "library(pacman);p_install('ggplot2')"
RUN R -e "library(pacman);p_install('rngtools')"
RUN R -e "library(pacman);p_install('NMF')"

RUN R -e "library(pacman);p_install('IRanges')"
RUN R -e "library(pacman);p_install('S4Vectors')"
RUN R -e "library(pacman);p_install('annotationDbi')"
RUN R -e "library(pacman);p_install('org.Hs.eg.db')"
Run R -e "library(pacman);p_install('rhdf5')"
RUN R -e "devtools::install_github('cmap/cmapR')"
RUN R -e "library(pacman);p_install('optparse')"
RUN R -e "library(pacman);p_install('verification')"
RUN R -e "library(pacman);p_install('gtools')"
RUN R -e "library(pacman);p_install('foreach')"
RUN R -e "library(pacman);p_install('doParallel')"
RUN R -e "library(pacman);p_install('magrittr')"

RUN R -e "library(pacman);p_install('maptools')"
RUN R -e "library(pacman);p_install('pheatmap')"
RUN R -e "library(pacman);p_install('WriteXLS')"
RUN R -e "library(pacman);p_install('RColorBrewer')"

RUN R -e "library(pacman);p_install('glue')"
RUN R -e "library(pacman);p_install('tidyselect')"
RUN R -e "library(pacman);p_install('dplyr')"

RUN R -e "library(pacman);p_install('Rtsne')"
RUN R -e "library(pacman);p_install('plotly')"
RUN R -e "library(pacman);p_install('rmarkdown')"
RUN R -e "library(pacman);p_install('vioplot')"
