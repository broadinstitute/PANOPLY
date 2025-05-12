FROM broadcptacdev/panoply_clumps_ptm:latest
MAINTAINER proteogenomics@broadinstitute.org

RUN apt-get update

COPY src /prot/proteomics/Projects/PGDAC/src





#### Install PyMol

RUN apt-get install -y pymol





#### Install PyMol
# Debian/Ubuntu/Mint
RUN apt-get install -y libglew-dev libpng-dev libfreetype6-dev libxml2-dev \
  libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev \
  libnetcdf-c++4-dev
RUN pip install cmake

# install from GitHub
RUN git clone https://github.com/schrodinger/pymol-open-source.git
RUN git clone https://github.com/rcsb/mmtf-cpp.git
RUN mv mmtf-cpp/include/mmtf* pymol-open-source/include/

# build PyMol
ENV prefix $HOME/pymol-open-source-build
ENV modules $prefix/modules
RUN cd pymol-open-source; python3 setup.py build install \
	--home=$prefix --install-lib=$modules

#add the module to python
RUN echo $modules > /usr/lib/python3/dist-packages/pymol.pth



#### Install PyMol

# Download Dependencies

RUN apt-get install -y libnetcdf-c++4-dev libglew-dev libglm-dev
RUN pip install cmake netCDF4

git clone https://github.com/schrodinger/pymol-open-source.git
git clone https://github.com/rcsb/mmtf-cpp.git
mv mmtf-cpp/include/mmtf* pymol-open-source/include/

prefix=$HOME/pymol-open-source-build

# Example for dependencies in non-standard places
# export PREFIX_PATH="$HOME/extra/glew-2.0.0:$HOME/extra/libpng-1.6.5:/opt/local"

python3 setup.py build install \
    --home=pymol-open-source/