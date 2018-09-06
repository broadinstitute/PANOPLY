FROM broadcptac/pgdac_cons_clust:1

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker for single-omics kmeans consensus clustering."


ENV wd /home/pgdac/
COPY . $wd
WORKDIR $wd

RUN rm -R /prot/proteomics/Projects/PGDAC/src/*
COPY src/* /prot/proteomics/Projects/PGDAC/src/

## samr not avauilable on CRAn since July 20, 2018
RUN R -e "install.packages('samr_2.0.tar.gz', repos=NULL, type='source', dependencies=T)"
RUN R -e "library(pacman);p_install('rjson')"
RUN R -e "library(pacman);p_install('GlobalOptions')"
RUN R -e "library(pacman);p_install('GetoptLong')"
RUN R -e "library(pacman);p_install('ComplexHeatmap')"
