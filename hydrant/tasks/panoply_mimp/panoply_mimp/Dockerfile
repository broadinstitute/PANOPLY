FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN apt-get update --allow-releaseinfo-change
RUN echo Y | apt-get install libxml2 libxml2-dev

RUN echo 'pacman::p_load (devtools, ComplexHeatmap, Biostrings)' \
        > /tmp/packages-devtools.R \
    && Rscript /tmp/packages-devtools.R

RUN echo 'if (! require (mimp)) { \
        devtools::install_github("omarwagih/rmimp", INSTALL_opts="--no-staged-install")};' \
      > /tmp/mimp-packages.R \
   && Rscript /tmp/mimp-packages.R

RUN R -e "install.packages(c('yaml', 'seqinr'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')"

RUN R -e "BiocManager::install('org.Hs.eg.db', update=FALSE)"
RUN R -e "BiocManager::install('EnsDb.Hsapiens.v79', update=FALSE)"

COPY src /prot/proteomics/Projects/PGDAC/src
