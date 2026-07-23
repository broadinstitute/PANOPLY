FROM proteomics/cosmo@sha256:581c53390e5b14d4e3e93185958ed358cbd27e03de575abe3542cb752dfeacad
MAINTAINER proteogenomics@broadinstitute.org

## set to stable repository
RUN sed -Ei 's/CRAN = .+?\), download/CRAN = \"packagemanager.rstudio.com\/cran\/2021-12-31\"\), download/' /usr/local/lib/R/etc/Rprofile.site # for rbase docker

## This base image's Rprofile.site adds BiocManager to `defaultPackages` (twice, redundantly),
## so it is auto-attached (library()'d) on every R startup, including every PSOCK cluster
## worker -- BiocManager's .onAttach() hook is what fetches bioconductor.org/config.yaml. This
## isn't gated by any BiocManager option, so remove the auto-attach itself; BiocManager is
## still installed and usable via BiocManager:: calls, it's just no longer force-loaded at
## startup.
RUN sed -i "/defaultPackages.*BiocManager/d" /usr/local/lib/R/etc/Rprofile.site
## Belt-and-suspenders: keep any other startup network call fast rather than a 60s hang.
RUN sed -i '1i options(timeout = 5)' /usr/local/lib/R/etc/Rprofile.site

COPY --from=broadcptacdev/panoply_common:latest /prot/proteomics/Projects/PGDAC/src /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src/cosmo

RUN R -e "BiocManager::install(c('cmapR', 'ComplexHeatmap'))"

RUN R -e "install.packages(c('pacman', 'optparse', 'rmarkdown', 'kableExtra'))"
