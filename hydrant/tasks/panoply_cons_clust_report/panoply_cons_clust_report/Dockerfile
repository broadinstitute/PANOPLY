FROM broadcptacdev/panoply_rmd_report:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN R -e "install.packages(c('yaml', 'stringr', 'tibble', 'tidyr', 'DT'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')"

COPY src /prot/proteomics/Projects/PGDAC/src
