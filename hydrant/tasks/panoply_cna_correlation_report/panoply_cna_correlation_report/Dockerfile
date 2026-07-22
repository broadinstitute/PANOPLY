FROM broadcptacdev/panoply_rmd_report:latest

LABEL maintainer="proteogenomics@broadinstitute.org"
LABEL description="Docker for generating Rmarkdown reports for different PANOPLY's CNA correlation module."

ENV wd /home/pgdac/

# Bundle jquery referenced by rmd-cna-analysis.r at build time,
# so the report renders offline (no CDN calls at runtime).
RUN mkdir -p $wd/lib \
 && curl -sSL -o $wd/lib/jquery-3.3.1.min.js https://code.jquery.com/jquery-3.3.1.min.js

COPY . $wd
WORKDIR $wd
