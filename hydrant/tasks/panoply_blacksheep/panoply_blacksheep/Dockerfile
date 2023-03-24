FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN R -e "BiocManager::install('blacksheepr', version='3.10')"

COPY src /prot/proteomics/Projects/PGDAC/src
