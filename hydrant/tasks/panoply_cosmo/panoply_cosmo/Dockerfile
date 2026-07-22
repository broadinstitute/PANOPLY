FROM proteomics/cosmo@sha256:581c53390e5b14d4e3e93185958ed358cbd27e03de575abe3542cb752dfeacad
MAINTAINER proteogenomics@broadinstitute.org

## set to stable repository
RUN sed -Ei 's/CRAN = .+?\), download/CRAN = \"packagemanager.rstudio.com\/cran\/2021-12-31\"\), download/' /usr/local/lib/R/etc/Rprofile.site # for rbase docker

## disable BiocManager's online version-diagnostics check (fetches bioconductor.org/config.yaml
## on every R startup, including each PSOCK cluster worker), which stalls ~60s per session
## when there is no internet access at runtime
ENV R_BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSTICS=FALSE
RUN echo 'options(BiocManager.online_version_diagnostics = FALSE)' >> /usr/local/lib/R/etc/Rprofile.site

COPY --from=broadcptacdev/panoply_common:latest /prot/proteomics/Projects/PGDAC/src /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src/cosmo

RUN R -e "BiocManager::install(c('cmapR', 'ComplexHeatmap'))"

RUN R -e "install.packages(c('pacman', 'optparse', 'rmarkdown', 'kableExtra'))"
