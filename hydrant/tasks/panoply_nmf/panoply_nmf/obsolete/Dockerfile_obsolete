FROM broadcptacdev/pgdac_mo_nmf:15

LABEL maintainer="proteogenomics@broadinstitute.org"
LABEL description="Docker for PANOPLY's multi-omic NMF clustering (v2.0)."

RUN R -e "pacman::p_install('pals')"

ENV wd /home/pgdac/
COPY . $wd
