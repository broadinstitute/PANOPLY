FROM broadcptacdev/panoply_clumps_ptm:latest
MAINTAINER proteogenomics@broadinstitute.org

# Install BLAST
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz
RUN tar -xf ncbi-blast-2.16.0+-x64-linux.tar.gz; rm ncbi-blast-2.16.0+-x64-linux.tar.gz
ENV PATH="$PATH:/ncbi-blast-2.16.0+/bin"


COPY src /prot/proteomics/Projects/PGDAC/src


