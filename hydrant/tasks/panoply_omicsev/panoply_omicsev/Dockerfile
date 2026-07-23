FROM proteomics/omicsev@sha256:235a707bc82551b1846488314f5b7e43127e19bd6135c2159981c76bac945fff
MAINTAINER proteogenomics@broadinstitute.org

## Removing BiocManager from `defaultPackages` to prevent loading in PSOCK cluster workers
## BiocManager's .onAttach() hook attempts to fetch bioconductor.org/config.yaml
## This causes a timeout error when run offline (e.g. on Manifold)
RUN sed -i "/defaultPackages.*BiocManager/d" /usr/local/lib/R/etc/Rprofile.site

COPY --from=broadcptacdev/panoply_common:latest /prot/proteomics/Projects/PGDAC/src /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src/omicsev

RUN R -e "BiocManager::install('cmapR')"

RUN R -e "install.packages(c('pacman', 'optparse'))"