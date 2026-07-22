FROM proteomics/omicsev@sha256:235a707bc82551b1846488314f5b7e43127e19bd6135c2159981c76bac945fff
MAINTAINER proteogenomics@broadinstitute.org

## disable BiocManager's online version-diagnostics check (fetches bioconductor.org/config.yaml
## on every R startup, including each PSOCK cluster worker), which stalls ~60s per session
## when there is no internet access at runtime
ENV R_BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSTICS=FALSE
RUN echo 'options(BiocManager.online_version_diagnostics = FALSE)' >> /usr/local/lib/R/etc/Rprofile.site

COPY --from=broadcptacdev/panoply_common:latest /prot/proteomics/Projects/PGDAC/src /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src/omicsev

RUN R -e "BiocManager::install('cmapR')"

RUN R -e "install.packages(c('pacman', 'optparse'))"