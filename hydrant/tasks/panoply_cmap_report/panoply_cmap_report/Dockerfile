FROM broadcptacdev/panoply_ssgsea_report
MAINTAINER proteogenomics@broadinstitute.org																														

## use same libdir as panoply_ssgsea_report
COPY src /home/pgdac/src/


RUN apt-get update
RUN apt-get install -y libpoppler-cpp-dev
RUN R -e "install.packages('pdftools')"



