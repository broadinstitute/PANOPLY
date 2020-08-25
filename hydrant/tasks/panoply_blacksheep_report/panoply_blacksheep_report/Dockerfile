FROM broadcptacdev/panoply_rmd_report:latest
MAINTAINER karen@broadinstitute.org

RUN R -e "install.packages('dplyr',\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')"

COPY src src
