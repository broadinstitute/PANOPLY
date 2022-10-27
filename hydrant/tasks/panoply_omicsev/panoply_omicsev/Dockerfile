FROM proteomics/omicsev@sha256:235a707bc82551b1846488314f5b7e43127e19bd6135c2159981c76bac945fff
MAINTAINER proteogenomics@broadinstitute.org

RUN mkdir /src

COPY src /src

RUN R -e "BiocManager::install('cmapR')"
