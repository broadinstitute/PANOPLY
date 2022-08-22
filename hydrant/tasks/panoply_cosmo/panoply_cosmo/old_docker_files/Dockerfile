FROM proteomics/cosmo:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN mkdir /opt/src

COPY src /opt/src

RUN R -e "BiocManager::install('cmapR')"
