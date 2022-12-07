FROM broadcptacdev/panoply_rmd_report:latest
MAINTAINER proteogenomics@broadinstitute.org

# the following packages are already installed in panoply_rmd_report
# RUN R -e "install.packages(c('yaml', 'stringr', 'tibble', 'tidyr', 'DT'),\
#     dependencies=TRUE,\
#     repos = 'http://cran.rstudio.com/')"

# needed for sankeyNetwork()
RUN R -e "install.packages(c('networkD3'),\
     dependencies=TRUE,\
     repos = 'http://cran.rstudio.com/')"																															

COPY src /prot/proteomics/Projects/PGDAC/src
