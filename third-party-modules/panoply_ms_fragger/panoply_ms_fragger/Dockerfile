FROM openjdk:8-jre-slim

RUN set -ex \
  && apt-get update && apt-get install -y --no-install-recommends zip \
  && rm -rf /var/lib/apt/lists/*

ADD src/MSFragger-2.1.jar /
COPY src/package.sh /usr/local/bin/

WORKDIR /MSFragger-2.1
