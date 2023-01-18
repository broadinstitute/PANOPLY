FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN R -e "install.packages('readr')"

# needed for sankeyNetwork()
RUN R -e "devtools::install_github(c('Displayr/flipU', 'Displayr/flipFormat', 'Displayr/flipTime'))"
# earlier version for proper formatting in .{tabset} HTML
RUN R -e "devtools::install_version('htmlwidgets', version = '1.5.4')"
RUN R -e "install.packages('/prot/proteomics/Projects/PGDAC/packages/flipTransformations_1.14.0.tar.gz')"																																																														

COPY src /prot/proteomics/Projects/PGDAC/src


