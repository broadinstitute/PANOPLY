FROM broadcptacdev/panoply_common:latest
MAINTAINER karen@broadinstitute.org

RUN R -e "install.packages(c('dplyr', 'stringr', 'tibble', 'yaml', 'BiocManager'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')" \
  && R -e "BiocManager::install('blacksheepr', version='3.10')"

COPY src /prot/proteomics/Projects/PGDAC/src
