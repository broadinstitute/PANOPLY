FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN set -ex \
    && cd /usr/local/bin \
    && apt-get update && apt-get install -y --no-install-recommends zip && apt-get install -y default-jre && apt install -y mono-devel && apt-get install -y r-base r-base-dev && apt install -y curl wget \
    && rm -rf /var/lib/apt/lists/* \
    && cd - \
    && cd /bin && ln -sf bash sh && cd - \
    && R -e "install.packages(c('optparse', 'pacman', 'yaml'), dependencies = TRUE, repos = 'http://cran.rstudio.com/')" 

RUN wget -qO /usr/local/bin/yq https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64
RUN chmod a+x /usr/local/bin/yq

RUN set -ex \
    && apt-get update \
    && wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh \
    && /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=/opt/conda/bin:$PATH
SHELL ["/bin/bash", "-c"]
RUN conda create --name fragpipe -y python=3.9
RUN echo "source activate fragpipe" > ~/.bashrc
ENV PATH /opt/conda/envs/fragpipe/bin:$PATH
RUN pip install Cython
RUN pip install matplotlib
RUN pip install easypqp
RUN pip install lxml

RUN set -ex \
    && mkdir ThermoRawFileParser && cd ThermoRawFileParser \
    && curl -LOk https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.1/ThermoRawFileParser1.4.1.zip \
    && unzip ThermoRawFileParser1.4.1.zip \
    && rm ThermoRawFileParser1.4.1.zip && cd -

ADD src/FragPipe-19.1.zip /
RUN unzip FragPipe-19.1.zip \
    && rm FragPipe-19.1.zip \
    && chmod 755 /fragpipe/bin/fragpipe

RUN set -ex \
    && cd /usr/local/bin \
    && curl -LOk https://github.com/Nesvilab/philosopher/releases/download/v4.8.0/philosopher_v4.8.0_linux_amd64.zip \
    && unzip philosopher_v4.8.0_linux_amd64.zip && rm philosopher_v4.8.0_linux_amd64.zip \
    && chmod 755 philosopher \
    && cd -

ADD src/IonQuant-1.8.10.zip /
RUN unzip IonQuant-1.8.10.zip && rm IonQuant-1.8.10.zip && chmod 755 /IonQuant-1.8.10/IonQuant-1.8.10.jar
 
ADD src/MSFragger-3.7.zip /
RUN unzip MSFragger-3.7.zip && rm MSFragger-3.7.zip && chmod 755 /MSFragger-3.7/MSFragger-3.7.jar

COPY src/collect_psm.sh /usr/local/bin/
COPY src/package.sh /usr/local/bin/
COPY src/mv.sh /etc/profile.d/
COPY src/get_fp_manifest.py /usr/local/bin/
