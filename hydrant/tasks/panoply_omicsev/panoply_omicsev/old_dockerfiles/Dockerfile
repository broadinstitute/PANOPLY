FROM proteomics/omicsev:latest
MAINTAINER svartany@broadinstitute.org

RUN mkdir /src

COPY src /src

RUN R -e "BiocManager::install('cmapR')"
