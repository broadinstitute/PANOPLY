FROM broadcptacdev/panoply_rmd_report:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN R -e "install.packages(c('yaml', 'tibble', 'tidyr'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')"

RUN wget https://raw.githubusercontent.com/broadinstitute/PANOPLY/dev/src/panoply_ssgsea_report/rmd-ssgsea-functions.R && mv rmd-ssgsea-functions.R /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src