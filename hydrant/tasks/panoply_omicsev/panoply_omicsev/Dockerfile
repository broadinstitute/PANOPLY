FROM proteomics/omicsev@sha256:235a707bc82551b1846488314f5b7e43127e19bd6135c2159981c76bac945fff
MAINTAINER proteogenomics@broadinstitute.org

## disable BiocManager's startup repository/version check, which stalls ~60s per R session
## (including each PSOCK cluster worker) when there is no internet access at runtime
RUN echo 'options(BiocManager.check_repositories = FALSE)' >> /usr/local/lib/R/etc/Rprofile.site

COPY --from=broadcptacdev/panoply_common:latest /prot/proteomics/Projects/PGDAC/src /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src/omicsev

RUN R -e "BiocManager::install('cmapR')"

RUN R -e "install.packages(c('pacman', 'optparse'))"