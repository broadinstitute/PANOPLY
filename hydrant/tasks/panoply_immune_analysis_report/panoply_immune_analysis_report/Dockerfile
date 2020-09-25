FROM broadcptacdev/panoply_common:latest
MAINTAINER karen@broadinstitute.org

RUN apt-get update
RUN apt-get -t unstable install -y pandoc

RUN echo 'pacman::p_load (devtools, ComplexHeatmap)' \
        > /tmp/packages-devtools.R \
    && Rscript /tmp/packages-devtools.R

RUN R -e "install.packages(c('tibble', 'dplyr', 'tidyr', 'rmarkdown', 'ggplot2', 'plotly'),\
    dependencies=TRUE,\
    repos = 'http://cran.rstudio.com/')"

COPY src /prot/proteomics/Projects/PGDAC/src
