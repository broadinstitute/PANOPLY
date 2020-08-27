FROM broadcptacdev/panoply_common:latest
MAINTAINER karen@broadinstitute.org

RUN R -e "install.packages(c('dplyr', 'rmarkdown'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')"

COPY src src
