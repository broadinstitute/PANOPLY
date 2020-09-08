FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:0.0.13

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker running gene-centric single-sample Gene Set Enrichment Analysis (ssGSEA) and PTM Signature Enrichment Analysis (PTM-SEA)".

ENV wd /home/pgdac/

COPY . $wd
WORKDIR $wd


## install R-packages
RUN R -e "install.packages('pacman')"
RUN R -e "library(pacman);p_install('rhdf5')"
#RUN R -e "devtools::install_github('cmap/cmapR')"
#RUN R -e "library(pacman);p_install('optparse')"
#RUN R -e "library(pacman);p_install('verification')"
#RUN R -e "library(pacman);p_install('gtools')"
#RUN R -e "library(pacman);p_install('foreach')"
#RUN R -e "library(pacman);p_install('doParallel')"
#RUN R -e "library(pacman);p_install('magrittr')"


#RUN rm src/d2*
#RUN rm src/phospho*
#RUN rm src/ssgsea-cli.R
