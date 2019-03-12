FROM broadcptac/pgdac_ssgsea:5

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker running gene-centric single-sample Gene Set Enrichment Analysis (ssGSEA) and PTM Signature Enrichment Analysis (PTM-SEA)".

ENV wd /home/pgdac/

COPY . $wd

WORKDIR $wd

RUN rm src/d2*
RUN rm src/phospho*
RUN rm src/ssgsea-cli.R
