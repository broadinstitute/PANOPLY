FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN R -e "install.packages(c('stringr', 'tibble', 'yaml', 'seqinr', 'remotes'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')" \
  && R -e "remotes::install_github('omarwagih/rmimp')"

COPY src /prot/proteomics/Projects/PGDAC/src
