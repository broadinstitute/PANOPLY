# Use an official Ubuntu as a parent image
FROM ubuntu:22.04

# Set noninteractive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

RUN set -ex; \
    cd /usr/local/bin; \
    apt-get update; \
    apt-get install -y --no-install-recommends zip default-jre mono-devel r-base r-base-dev curl wget; \
    rm -rf /var/lib/apt/lists/*; \
    cd -; \
    cd /bin && ln -sf bash sh && cd -; \
    R -e "install.packages(c('optparse', 'pacman', 'yaml'), dependencies = TRUE, repos = 'http://cran.rstudio.com/')"; \
    
    wget -qO /usr/local/bin/yq https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64; \
    chmod a+x /usr/local/bin/yq; \
    
    # miniconda
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh; \
    /bin/bash ~/miniconda.sh -b -p /opt/conda; \
    
    # python and fragpipe
    PATH=/opt/conda/bin:$PATH; \
    conda create --name fragpipe -y python=3.9; \
    echo "source activate fragpipe" > ~/.bashrc; \
    PATH=/opt/conda/envs/fragpipe/bin:$PATH; \
    pip install Cython matplotlib easypqp lxml; \
    
    # ThermoRawFileParser
    cd /; \
    mkdir ThermoRawFileParser && cd ThermoRawFileParser; \
    curl -LOk https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.3/ThermoRawFileParser1.4.3.zip; \
    unzip ThermoRawFileParser1.4.3.zip; \
    rm ThermoRawFileParser1.4.3.zip; \

    # FragPipe
    cd /; \
    curl -LOk https://github.com/Nesvilab/FragPipe/releases/download/20.0/FragPipe-20.0.zip; \
    unzip FragPipe-20.0.zip; \
    rm FragPipe-20.0.zip; \
    chmod 755 /fragpipe/bin/fragpipe; \
    
    # philosopher_v5
    cd /usr/local/bin; \
    curl -LOk https://github.com/Nesvilab/philosopher/releases/download/v5.0.0/philosopher_v5.0.0_linux_amd64.zip; \
    unzip philosopher_v5.0.0_linux_amd64.zip; \
    rm philosopher_v5.0.0_linux_amd64.zip; \
    chmod 755 philosopher

# IonQuant -- had to add here as SW not avalaible online for free
ADD src/IonQuant-1.9.8.zip /
RUN unzip IonQuant-1.9.8.zip && rm IonQuant-1.9.8.zip && chmod 755 /IonQuant-1.9.8/IonQuant-1.9.8.jar
 
# MSFragger -- had to add here as SW not avalaible online for free
ADD src/MSFragger-3.8.zip /
RUN unzip MSFragger-3.8.zip && rm MSFragger-3.8.zip && chmod 755 /MSFragger-3.8/MSFragger-3.8.jar

# copy files in the src directory to usr/local/bin
COPY src/collect_psm.sh src/package.sh src/mv.sh src/get_fp_manifest.py /usr/local/bin/

# Set the PATH environment variable
ENV PATH=/opt/conda/bin:/opt/conda/envs/fragpipe/bin:$PATH

# Set the default shell to /bin/bash
SHELL ["/bin/bash", "-c"]
