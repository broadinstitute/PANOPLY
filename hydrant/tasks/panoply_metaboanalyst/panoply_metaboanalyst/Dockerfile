FROM rocker/r-ver:4.1
MAINTAINER proteogenomics@broadinstitute.org

RUN apt-get update
RUN apt-get install -y libcairo2-dev libnetcdf-dev libxml2 libxt-dev libssl-dev libbz2-dev libglpk-dev curl texlive-latex-recommended
RUN R -e 'install.packages(c("optparse", "tidyverse", "rjson"))'
#RUN R -e 'devtools::install_github("metabolomicsworkbench/RefMet")'

COPY src /prot/proteomics/Projects/PGDAC/src

# install MetaboAnalyst dependencies & MetaboAnalyst
RUN Rscript /prot/proteomics/Projects/PGDAC/src/install_dependencies.R
RUN R -e 'BiocManager::install("cmapR")'
RUN R -e 'devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)'