FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN R -e "install.packages(c('yaml', 'BiocManager'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')" \
  && R -e "BiocManager::install('blacksheepr', version='3.10')"

COPY src /prot/proteomics/Projects/PGDAC/src
