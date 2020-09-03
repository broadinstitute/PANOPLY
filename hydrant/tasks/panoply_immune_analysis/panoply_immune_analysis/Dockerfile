FROM broadcptacdev/panoply_common:latest
MAINTAINER rkothadi@broadinstitute.org

RUN apt-get update 
RUN echo Y | apt-get install libxml2 libxml2-dev

RUN echo 'pacman::p_load (devtools, ComplexHeatmap)' \
        > /tmp/packages-devtools.R \
    && Rscript /tmp/packages-devtools.R

RUN echo 'if (! require (estimate)) { \
        install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)}; \
      if (! require (XML)) { \
        install.packages("XML", repos = "http://www.omegahat.net/R")}; \ 
      if (! require (xCell)) { \
        devtools::install_github("dviraran/xCell")}; \
      if (! require (ImmuneSubtypeClassifier)) { \
        devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")};' \
      > /tmp/immune-packages.R \
   && Rscript /tmp/immune-packages.R

        
COPY src /prot/proteomics/Projects/PGDAC/src
