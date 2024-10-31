FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN apt-get update
RUN apt-get install -y libbz2-dev

# install additional dependencies
RUN R -e "BiocManager::install(c('globaltest', 'ActivePathways', 'org.Hs.eg.db', 'KEGGgraph'))"
RUN R -e "install.packages(c('qs','statmod', 'igraph', 'SimDesign'))"

COPY src /prot/proteomics/Projects/PGDAC/src