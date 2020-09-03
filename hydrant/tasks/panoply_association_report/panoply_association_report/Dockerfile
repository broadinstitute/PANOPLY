FROM broadcptacdev/panoply_common:latest
MAINTAINER karen@broadinstitute.org

RUN apt-get update
RUN apt-get -t unstable install -y pandoc

RUN R -e "install.packages(c('yaml', 'tibble', 'dplyr', 'tidyr', 'rmarkdown', 'ggplot2', 'plotly', 'glue'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')"

RUN wget https://raw.githubusercontent.com/broadinstitute/PANOPLY/issue-add-mo-nmf/src/panoply_ssgsea_report/rmd-ssgsea-functions.R && mv rmd-ssgsea-functions.R /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src