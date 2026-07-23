FROM proteomics/cosmo@sha256:581c53390e5b14d4e3e93185958ed358cbd27e03de575abe3542cb752dfeacad
MAINTAINER proteogenomics@broadinstitute.org

## set to stable repository
RUN sed -Ei 's/CRAN = .+?\), download/CRAN = \"packagemanager.rstudio.com\/cran\/2021-12-31\"\), download/' /usr/local/lib/R/etc/Rprofile.site # for rbase docker

## Removing BiocManager from `defaultPackages` to prevent loading in PSOCK cluster workers
## BiocManager's .onAttach() hook attempts to fetch bioconductor.org/config.yaml
## This causes a timeout error when run offline (e.g. on Manifold)
RUN sed -i "/defaultPackages.*BiocManager/d" /usr/local/lib/R/etc/Rprofile.site

COPY --from=broadcptacdev/panoply_common:latest /prot/proteomics/Projects/PGDAC/src /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src/cosmo

RUN R -e "BiocManager::install(c('cmapR', 'ComplexHeatmap'))"

RUN R -e "install.packages(c('pacman', 'optparse', 'rmarkdown', 'kableExtra'))"
