#FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:0.0.13
FROM broadcptacdev/panoply_common:74e21ac

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker running gene-centric single-sample Gene Set Enrichment Analysis (ssGSEA) and PTM Signature Enrichment Analysis (PTM-SEA)".

#ENV wd /home/pgdac/
ENV wd /home/jupyter-user


COPY . $wd
WORKDIR $wd

RUN ls -al
#RUN sudo chown jupyter-user:users src
RUN ls src

# get ssGSEA scripts from GitHub                                                
#RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/ssGSEA2.0.R
#RUN dos2unix ssGSEA2.0.R
#RUN cp ssGSEA2.0.R src

#RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/parse_yaml.R
#RUN dos2unix parse_yaml.R

#RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/ssgsea-cli.R
#RUN dos2unix ssgsea-cli.R



## install R-packages
RUN R -e "install.packages('pacman')"
#RUN R -e "pacman::p_install('rhdf5')"
#RUN R -e "devtools::install_github('cmap/cmapR')"
#RUN R -e "pacman::p_install('optparse')"
#RUN R -e "pacman::p_install('verification')"
#RUN R -e "pacman::p_install('gtools')"
#RUN R -e "pacman::p_install('foreach')"
#RUN R -e "pacman::p_install('doParallel')"
#RUN R -e "pacman::p_install('magrittr')"
