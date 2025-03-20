FROM python:3.9
MAINTAINER proteogenomics@broadinstitute.org

RUN apt-get update
RUN apt-get install -y pv zstd pigz

RUN pip config set global.index-url https://packagemanager.posit.co/pypi/2025-01-01/simple
RUN pip config set global.trusted-host packagemanager.posit.co


# # Install Clumps-PTM
# RUN pip3 install clumps-ptm


# Install Clumps-PTM Dependencies
RUN pip3 install prody matplotlib pyarrow # buildkit

# Install CLUMPS-PTM from Git Repo
COPY src/CLUMPS-PTM/ CLUMPS-PTM/
RUN cd CLUMPS-PTM/; python setup.py install


# Install Other Packages
RUN pip3 install plotting requests cmapPy pyyaml

