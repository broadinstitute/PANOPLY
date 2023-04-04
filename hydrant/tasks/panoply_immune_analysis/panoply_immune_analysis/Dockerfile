FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN apt-get update 
RUN echo Y | apt-get install libxml2 libxml2-dev

RUN echo 'if (! require (xCell)) { \
        devtools::install_github("dviraran/xCell")}; \
      if (! require (ImmuneSubtypeClassifier)) { \
        devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")};' \
      > /tmp/immune-packages.R \
   && Rscript /tmp/immune-packages.R

        
COPY src /prot/proteomics/Projects/PGDAC/src
