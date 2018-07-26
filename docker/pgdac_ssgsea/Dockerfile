FROM broadcptac/r-util:1

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker running gene-centric single-sample Gene Set Enrichment Analysis (ssGSEA) and site-specific PTM Signature Enrichment Analysis (ssPSEA)".

ENV wd /home/pgdac/

COPY . $wd
#COPY src src/

WORKDIR $wd

## install R-packages
RUN R -e "install.packages('pacman')"
Run R -e "library(pacman);p_install('rhdf5')"
RUN R -e "devtools::install_github('cmap/cmapR')"
RUN R -e "library(pacman);p_install('optparse')"
RUN R -e "library(pacman);p_install('verification')"
RUN R -e "library(pacman);p_install('gtools')"
RUN R -e "library(pacman);p_install('foreach')"
RUN R -e "library(pacman);p_install('doParallel')"
RUN R -e "library(pacman);p_install('magrittr')"


