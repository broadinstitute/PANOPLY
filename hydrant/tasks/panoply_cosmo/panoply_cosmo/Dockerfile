FROM proteomics/cosmo@sha256:581c53390e5b14d4e3e93185958ed358cbd27e03de575abe3542cb752dfeacad
MAINTAINER proteogenomics@broadinstitute.org

## set to stable repository
RUN sed -Ei 's/CRAN = .+?\), download/CRAN = \"packagemanager.rstudio.com\/cran\/2021-12-31\"\), download/' /usr/local/lib/R/etc/Rprofile.site # for rbase docker

## Whatever in this base image's Rprofile.site triggers a bioconductor.org/config.yaml fetch on
## every R startup (including each PSOCK cluster worker) runs before any options we *append* to
## the file, so appending BiocManager-specific switches (check_repositories,
## online_version_diagnostics) has proven ineffective -- the call already ran by the time our
## line is reached. Instead, PREPEND a short, generic base-R connection timeout so that whatever
## network call fires first fails in a few seconds instead of hanging the default 60s per R
## session -- this works regardless of which package/mechanism makes the call.
RUN sed -i '1i options(timeout = 5)' /usr/local/lib/R/etc/Rprofile.site
ENV R_BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSTICS=FALSE
RUN echo 'options(BiocManager.online_version_diagnostics = FALSE)' >> /usr/local/lib/R/etc/Rprofile.site

COPY --from=broadcptacdev/panoply_common:latest /prot/proteomics/Projects/PGDAC/src /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src/cosmo

RUN R -e "BiocManager::install(c('cmapR', 'ComplexHeatmap'))"

RUN R -e "install.packages(c('pacman', 'optparse', 'rmarkdown', 'kableExtra'))"
