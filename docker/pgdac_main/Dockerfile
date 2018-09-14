FROM broadcptac/r-util:2
MAINTAINER manidr@broadinstitute.org

COPY src /prot/proteomics/Projects/PGDAC/src
COPY data /prot/proteomics/Projects/PGDAC/data
COPY R-utilities /prot/proteomics/Projects/R-utilities

## install additional libraries required for PGDAC pipeline components
RUN R -e "library(pacman);p_install('optparse')"
RUN R -e "library(pacman);p_install('NMF')"
RUN R -e "library(pacman);p_install('factoextra')"
RUN R -e "library(pacman);p_install('glue')"
RUN R -e "library(pacman);p_install('magrittr')"
Run R -e "library(pacman);p_install('rhdf5')"
RUN R -e "devtools::install_github('cmap/cmapR')"
RUN R -e "library(pacman);p_install('doParallel')"
RUN R -e "library(pacman);p_install('foreach')"
RUN R -e "library(pacman);p_install('pheatmap')"
RUN R -e "library(pacman);p_install('RColorBrewer')"
RUN R -e "library(pacman);p_install('NbClust')"
# samr not avauilable on CRAN since July 20, 2018
RUN R -e "install.packages('packages/samr_2.0.tar.gz', repos=NULL, type='source', dependencies=T)"
RUN R -e "library(pacman);p_install('rjson')"
RUN R -e "library(pacman);p_install('GlobalOptions')"
RUN R -e "library(pacman);p_install('GetoptLong')"
RUN R -e "library(pacman);p_install('ComplexHeatmap')"
RUN R -e "library(pacman);p_install('dplyr')"
