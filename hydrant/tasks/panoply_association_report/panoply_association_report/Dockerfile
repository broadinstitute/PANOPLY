FROM broadcptacdev/panoply_common:latest
MAINTAINER karen@broadinstitute.org

RUN R -e "install.packages(c('dplyr', 'rmarkdown'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')"

RUN wget https://raw.githubusercontent.com/broadinstitute/PANOPLY/issue-add-mo-nmf/src/panoply_ssgsea_report/rmd-ssgsea-functions.R && mv rmd-ssgsea-functions.R /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src