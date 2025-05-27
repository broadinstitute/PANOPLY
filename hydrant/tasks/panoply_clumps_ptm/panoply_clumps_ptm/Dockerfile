FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org
CMD ["/usr/bin/python3"]

RUN apt-get update
RUN apt-get install -y parallel


# Install PyMol (and python 3.7.3) and Pip / Wheel
RUN apt-get install -y pymol python3-pip
RUN ln -s /usr/bin/python3 /usr/bin/python # add python pointer to python3


# set posit package management to a stable repository
RUN pip3 config set global.index-url https://packagemanager.posit.co/pypi/2025-01-01/simple
RUN pip3 config set global.trusted-host packagemanager.posit.co

RUN pip3 install --upgrade pip
RUN pip install setuptools wheel build tomli




# Install CLUMPS-PTM dependencies
RUN apt-get install -y python3-statsmodels
RUN pip install packaging "Cython>=3.0" pyarrow

# Install CLUMPS-PTM from Git Repo
ADD https://github.com/getzlab/CLUMPS-PTM.git /CLUMPS-PTM
RUN pip install -e CLUMPS-PTM/


# Install Other Python Packages
RUN pip install plotting requests cmapPy pyyaml pprint3x "numexpr>=2.7.0"





COPY src /prot/proteomics/Projects/PGDAC/src