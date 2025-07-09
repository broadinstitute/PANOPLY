FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN set -ex \
    && apt-get update \
    && apt-get -y install dotnet-sdk-8.0 \
    && apt-get install -y --no-install-recommends zip \
    && apt-get install -y clang \
    && rm -rf /var/lib/apt/lists/*
# NB: Ubuntu does not include zip
#     LightGBM (used in Spectronaut) requires clang

COPY packages/Spectronaut_20.1.250624.92449.deb /packages/Spectronaut_20.1.250624.92449.deb
RUN dpkg -i /packages/Spectronaut_20.1.250624.92449.deb

RUN spectronaut -activate #license key