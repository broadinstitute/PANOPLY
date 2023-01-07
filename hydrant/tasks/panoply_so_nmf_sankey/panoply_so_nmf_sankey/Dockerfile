FROM broadcptacdev/panoply_rmd_report:latest
MAINTAINER proteogenomics@broadinstitute.org

# the following packages are already installed in panoply_rmd_report
# RUN R -e "install.packages(c('yaml', 'stringr', 'tibble', 'tidyr', 'DT'),\
#     dependencies=TRUE,\
#     repos = 'http://cran.rstudio.com/')"

# needed for sankeyNetwork()
RUN R -e "require(devtools);\
          install_github(c('Displayr/flipU', 'Displayr/flipFormat', 'Displayr/flipTime'))"
# earlier version for proper formatting in tabset
RUN R -e "require(devtools);\
          install_version('htmlwidgets', version = '1.5.4')"
RUN R -e "install.packages('/prot/proteomics/Projects/PGDAC/packages/flipTransformations_1.14.0.tar.gz')"																																																														

COPY src /prot/proteomics/Projects/PGDAC/src


