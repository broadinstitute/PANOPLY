FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN set -ex \
    && cd /usr/local/bin \
    && apt-get update && apt install -y curl libgomp1 wget gnupg software-properties-common zip mono-devel \
    && rm -rf /var/lib/apt/lists/* \
    && curl -LOk https://github.com/vdemichev/DiaNN/releases/download/1.8.1/diann_1.8.1.deb \
    && dpkg -i ./diann_1.8.1.deb \
    && ln -sf /usr/diann/1.8.1/diann-1.8.1 /usr/bin/diann \   
    && cd - \
    && cd /bin && ln -sf bash sh && cd -

RUN set -ex \
    && mkdir ThermoRawFileParser && cd ThermoRawFileParser \
    && curl -LOk https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.0/ThermoRawFileParser.zip \
    && unzip ThermoRawFileParser.zip \
    && rm ThermoRawFileParser.zip && cd -


