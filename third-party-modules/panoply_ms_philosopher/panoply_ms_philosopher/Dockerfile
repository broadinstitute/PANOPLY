FROM broadgdac/firecloud-ubuntu:16.04

RUN set -ex \
    && cd /usr/local/bin \
    && apt-get update && apt-get install -y --no-install-recommends zip \
    && curl -LOk https://github.com/Nesvilab/philosopher/releases/download/v2.0.0/philosopher_v2.0.0_linux_amd64.zip \
    && unzip philosopher_v2.0.0_linux_amd64.zip \
    && chmod 755 philosopher \
    && cd - \
    && cd /bin && ln -sf bash sh && cd -

COPY src/mv.sh /etc/profile.d/

CMD ["philosopher"]
