FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

# needed for sankeyNetwork()
RUN R -e "library(pacman);p_install('networkD3')"
RUN R -e "devtools::install_github(c('Displayr/flipU', 'Displayr/flipFormat', 'Displayr/flipTime'))"
RUN R -e "install.packages('/prot/proteomics/Projects/PGDAC/packages/flipTransformations_1.14.0.tar.gz')"

# earlier version for proper formatting in .{tabset} HTML
RUN R -e "devtools::install_version('htmlwidgets', version = '1.5.4')"


COPY src /prot/proteomics/Projects/PGDAC/src


