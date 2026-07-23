FROM proteomics/omicsev@sha256:235a707bc82551b1846488314f5b7e43127e19bd6135c2159981c76bac945fff
MAINTAINER proteogenomics@broadinstitute.org

## This base image's Rprofile.site adds BiocManager to `defaultPackages`, so it is
## auto-attached (library()'d) on every R startup, including every PSOCK cluster worker --
## BiocManager's .onAttach() hook is what fetches bioconductor.org/config.yaml. This isn't
## gated by any BiocManager option, so remove the auto-attach itself; BiocManager is still
## installed and usable via BiocManager:: calls, it's just no longer force-loaded at startup.
RUN sed -i "/defaultPackages.*BiocManager/d" /usr/local/lib/R/etc/Rprofile.site

COPY --from=broadcptacdev/panoply_common:latest /prot/proteomics/Projects/PGDAC/src /prot/proteomics/Projects/PGDAC/src

COPY src /prot/proteomics/Projects/PGDAC/src/omicsev

RUN R -e "BiocManager::install('cmapR')"

RUN R -e "install.packages(c('pacman', 'optparse'))"

## Belt-and-suspenders: keep any *runtime* startup network call fast rather than a 60s hang.
## Applied only after the package installs above, since Rprofile.site is sourced by every R
## invocation including these build steps -- a 5s timeout here previously caused
## BiocManager::install() to silently drop packages (e.g. cmapR) that took >5s to download,
## without failing the build (install.packages() only warns on a per-package failure).
RUN sed -i '1i options(timeout = 5)' /usr/local/lib/R/etc/Rprofile.site