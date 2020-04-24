FROM broadcptac/pgdac_mo_nmf:11

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker for multi-omic NMF clustering."

ENV wd /home/pgdac/

RUN rm -R src

COPY . $wd

RUN R -e "pacman::p_install('ggpubr')"

