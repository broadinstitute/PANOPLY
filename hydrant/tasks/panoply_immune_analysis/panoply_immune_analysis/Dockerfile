FROM broadcptacdev/panoply_common:70df9e7
MAINTAINER rkothadi@broadinstitute.org

RUN apt-get update \
    && apt-get install libxml2 libxml2-dev

RUN echo 'pacman::p_load (devtools, ComplexHeatmap)' \
        > /tmp/packages-devtools.R \
    && Rscript /tmp/packages-devtools.R

RUN echo 'if (! require (estimate)) { \
        install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)}; \
      if (! require (xCell)) { \
        devtools::install_github("dviraran/xCell")}; \
      if (! require (ImmuneSubtypeClassifier)) { \
        devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")};' \
      > /tmp/immune-packages.R \
   && Rscript /tmp/immune-packages.R

        
COPY src /prot/proteomics/Projects/PGDAC/src
