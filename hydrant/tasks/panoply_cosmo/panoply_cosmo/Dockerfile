FROM proteomics/cosmo@sha256:581c53390e5b14d4e3e93185958ed358cbd27e03de575abe3542cb752dfeacad
MAINTAINER proteogenomics@broadinstitute.org

COPY --from=broadcptacdev/panoply_common:latest /prot/proteomics/Projects/PGDAC/src /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src/cosmo

RUN R -e "BiocManager::install(c('cmapR', 'ComplexHeatmap'))"

RUN R -e "install.packages(c('pacman', 'optparse', 'rmarkdown', 'kableExtra'))"
