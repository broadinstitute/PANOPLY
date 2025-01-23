FROM python:3.9
MAINTAINER proteogenomics@broadinstitute.org

RUN pip config set global.index-url https://packagemanager.posit.co/pypi/2025-01-01/simple
RUN pip config set global.trusted-host packagemanager.posit.co

RUN pip3 install clumps-ptm


COPY src /prot/proteomics/Projects/PGDAC/src


