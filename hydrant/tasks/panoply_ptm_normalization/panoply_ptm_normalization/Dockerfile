FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

## change CRAN repository to April 2023, to allow use of tidyr::separate_longer_delim() 
# RUN sed -i 's/2020-12-31/2023-04-03/' /etc/R/Rprofile.site # for rbase docker
RUN sed -i 's/2020-12-31/2023-04-03/' /usr/local/lib/R/etc/Rprofile.site # for rocker docker
## update tidyverse package, to allow use of tidyr::separate_longer_delim() 
RUN R -e 'install.packages( "tidyverse" )'

COPY src /prot/proteomics/Projects/PGDAC/src
