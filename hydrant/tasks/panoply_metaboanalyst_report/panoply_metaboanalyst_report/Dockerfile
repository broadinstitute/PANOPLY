FROM broadcptacdev/panoply_rmd_report:latest

LABEL maintainer="proteogenomics@broadinstitute.org"
LABEL description="Docker for generating Rmarkdown reports for different PGDAC modules."

ENV wd /prot/proteomics/Projects/PGDAC/src

# Bundle JS/CSS assets referenced by metaboanalyst_rmd.rmd at build time,
# so the report renders offline (no CDN calls at runtime).
RUN mkdir -p $wd/lib \
 && curl -sSL -o $wd/lib/jquery.min.js https://ajax.googleapis.com/ajax/libs/jquery/2.1.3/jquery.min.js \
 && curl -sSL -o $wd/lib/d3.v4.js https://d3js.org/d3.v4.js \
 && curl -sSL -o $wd/lib/dataTables.dataTables.css https://cdn.datatables.net/2.1.8/css/dataTables.dataTables.css \
 && curl -sSL -o $wd/lib/dataTables.js https://cdn.datatables.net/2.1.8/js/dataTables.js

COPY src /prot/proteomics/Projects/PGDAC/src
WORKDIR $wd