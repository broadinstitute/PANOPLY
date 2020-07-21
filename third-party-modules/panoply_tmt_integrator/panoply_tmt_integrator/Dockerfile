FROM openjdk:8-jre-slim

RUN set -ex \
  && apt-get update && apt-get install -y --no-install-recommends zip unzip \
  && rm -rf /var/lib/apt/lists/* \
  && mkdir -p TMTIntegrator_v1.1.0/results

ADD src/TMTIntegrator_v1.1.0.jar /
COPY src/package.sh /usr/local/bin/
COPY src/collect_psm.sh /usr/local/bin/

WORKDIR /TMTIntegrator_v1.1.0
