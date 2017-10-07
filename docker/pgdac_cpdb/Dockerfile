FROM broadcptac/r-util:1
#FROM bioconductor/devel_base2

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker for generating customized protein sequence database from genomics data using customProDB."

ENV wd /home/pgdac/

COPY . $wd
COPY src src/

WORKDIR $wd

RUN apt-get update
## required for R-pacakge 'XML'
RUN apt-get -t unstable install -y libxml2-dev
## rmarkdown
RUN apt-get -t unstable install -y pandoc


## install R-packages
RUN R -e "install.packages('pacman')"
RUN R -e "library(pacman);p_install('rmarkdown')"
RUN R -e "library(pacman);p_install('seqinr')"
RUN R -e "library(pacman);p_install('BSgenome.Hsapiens.UCSC.hg19')"
RUN R -e "library(pacman);p_install('WriteXLS')"
RUN R -e "library(pacman);p_install('shiny')"
RUN R -e "library(pacman);p_install('plotly')"
RUN R -e "library(pacman);p_install('customProDB')"
RUN R -e "library(pacman);p_install('GenomicFeatures')"
RUN R -e "library(pacman);p_install('AnnotationDbi')"
RUN R -e "install.packages('customProDBBI_1.13.0.tar.gz', repos=NULL, type='source', dependencies=T)"


## remove tar.gz of 'customProDBBI'
RUN rm customProDBBI_1.13.0.tar.gz
RUN rm -R src

