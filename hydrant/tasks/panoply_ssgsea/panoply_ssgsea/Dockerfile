#FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:0.0.13
FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

LABEL description="Docker running gene-centric single-sample Gene Set Enrichment Analysis (ssGSEA) and PTM Signature Enrichment Analysis (PTM-SEA)".

ENV wd /home/pgdac/
#ENV wd /home/jupyter-user
#RUN sudo chown jupyter-user:users src

COPY . $wd
WORKDIR $wd

RUN apt-get update \
	&& apt-get install -y wget \
	&& apt-get install -y dos2unix \
	&& apt-get install -y nano

RUN dos2unix src/*

# get ssGSEA scripts from GitHub                                                
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/ssGSEA2.0.R && dos2unix ssGSEA2.0.R && mv ssGSEA2.0.R src
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/parse_yaml_ssgsea.R && dos2unix parse_yaml_ssgsea.R && mv parse_yaml_ssgsea.R src
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/io.R && dos2unix io.R && mv io.R src
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/utils.R && dos2unix utils.R && mv utils.R src
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/ssgsea-cli.R && dos2unix ssgsea-cli.R && chmod +x ssgsea-cli.R

## install R-packages
#RUN R -e "install.packages('pacman')"
#RUN R -e "pacman::p_install('rhdf5')"
#RUN R -e "devtools::install_github('cmap/cmapR')"
#RUN R -e "pacman::p_install('optparse')"
#RUN R -e "pacman::p_install('verification')"
#RUN R -e "pacman::p_install('gtools')"
#RUN R -e "pacman::p_install('foreach')"
#RUN R -e "pacman::p_install('doParallel')"
#RUN R -e "pacman::p_install('magrittr')"
