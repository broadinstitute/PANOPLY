FROM bioconductor/devel_base2

LABEL maintainer="karsten@broadinstitute.org"
LABEL description="Docker for generating customized protein sequence database from genomics data using customProDB."

ENV wd /home/pgdac/

COPY . $wd
COPY src src/

WORKDIR $wd

##VOLUME /share/

## install R-packages
RUN R -e "rmarkdown"
RUN R -e "biocLite('seqinr')"
RUN R -e "biocLite('customProDB')"
RUN R -e "biocLite('BSgenome.Hsapiens.UCSC.hg19')"
RUN R -e "install.packages('customProDBBI_1.13.0.tar.gz', repos=NULL, type='source', dependencies=T)"
RUN R -e "install.packages('WriteXLS')"


## remove tar.gz of 'customProDBBI'
RUN rm customProDBBI_1.13.0.tar.gz

