FROM broadcptacdev/panoply_common:latest

LABEL maintainer="proteogenomics@broadinstitute.org"
LABEL description="Docker for generating Rmarkdown reports for different Panoply modules."

ENV wd /home/pgdac/

COPY . $wd
WORKDIR $wd

## install R-packages
RUN R -e "devtools::install_github('cmap/morpheus.R')"
RUN R -e "library(pacman);p_install('DT')"
RUN R -e "library(pacman);p_install('gtools')"
RUN R -e "library(pacman);p_install('rmarkdown')"
RUN R -e "library(pacman);p_install('WriteXLS')"
RUN R -e "library(pacman);p_install('shiny')"
RUN R -e "library(pacman);p_install('plotly')"
RUN R -e "library(pacman);p_install('seriation')"
RUN R -e "library(pacman);p_install('kableExtra')"


## install missing fonts and images required for morpheus
# RUN mv fonts/* /usr/local/lib/R/site-library/morpheus/htmlwidgets/lib/css/fonts/
# RUN mv images/* /usr/local/lib/R/site-library/morpheus/htmlwidgets/lib/css/images/
# RUN rm -R fonts
# RUN rm -R images
