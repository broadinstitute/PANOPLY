FROM python:3.9
MAINTAINER proteogenomics@broadinstitute.org

RUN pip config set global.index-url https://packagemanager.posit.co/pypi/2025-01-01/simple
RUN pip config set global.trusted-host packagemanager.posit.co

# Install Clumps-PTM
RUN pip3 install clumps-ptm

# Install Other Packages
RUN pip3 install plotting

# Install BLAST
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz
RUN tar -xf ncbi-blast-2.16.0+-x64-linux.tar.gz
ENV PATH="$PATH:/ncbi-blast-2.16.0+/bin"


COPY src /prot/proteomics/Projects/PGDAC/src


