FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN set -ex \
    && apt-get update \
    && apt-get -y install dotnet-sdk-6.0 \
    && apt-get install -y --no-install-recommends zip \
    && apt-get install -y clang \
    && rm -rf /var/lib/apt/lists/*
# NB: Ubuntu does not include zip
#     LightGBM (used in Spectronaut) requires clang

COPY packages/Spectronaut_19.0.240606.62635.deb /packages/Spectronaut_19.0.240606.62635.deb
RUN dpkg -i /packages/Spectronaut_19.0.240606.62635.deb
