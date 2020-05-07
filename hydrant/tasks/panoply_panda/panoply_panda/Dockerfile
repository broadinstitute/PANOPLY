FROM ubuntu:bionic
MAINTAINER rkothadi@broadinstitute.org

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update
RUN apt-get install -y curl python gcc python-dev python-setuptools libffi-dev
RUN apt-get install -y python-pip vim
RUN pip install gsutil
RUN pip install firecloud

COPY bin /panda/bin
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" \ 
      >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com \ 
      --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt update
RUN apt install -y r-base r-base-dev
