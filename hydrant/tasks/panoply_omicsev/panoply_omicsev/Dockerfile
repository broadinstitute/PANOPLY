FROM proteomics/omicsev@sha256:235a707bc82551b1846488314f5b7e43127e19bd6135c2159981c76bac945fff
MAINTAINER proteogenomics@broadinstitute.org

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

COPY src /prot/proteomics/Projects/PGDAC/src/omicsev

RUN R -e "BiocManager::install('cmapR')"

RUN R -e "install.packages(c('pacman', 'optparse'))"