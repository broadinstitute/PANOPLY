FROM broadcptac/pgdac_rmd:3

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker for generating Rmarkdown reports for different PANOPLY modules."

ENV wd /home/pgdac/

RUN R -e "library(pacman);p_install('yaml')"

COPY . $wd
WORKDIR $wd
