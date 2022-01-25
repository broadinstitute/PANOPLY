FROM broadgdac/firecloud-ubuntu:16.04

RUN set -ex \
    && cd /usr/local/bin \
    && apt-get update && apt-get install -y --no-install-recommends zip && apt-get install -y default-jre && apt install -y mono-devel \
    && rm -rf /var/lib/apt/lists/* \
    && curl -LOk https://github.com/Nesvilab/philosopher/releases/download/v4.1.1/philosopher_v4.1.1_linux_amd64.zip \
    && unzip philosopher_v4.1.1_linux_amd64.zip \
    && chmod 755 philosopher \
    && cd - \
    && cd /bin && ln -sf bash sh && cd -

ADD src/MSFragger-3.4.zip /
RUN unzip MSFragger-3.4.zip
RUN rm MSFragger-3.4.zip

RUN curl -LOk https://github.com/Nesvilab/TMT-Integrator/releases/download/3.1.0/TMTIntegrator-v3.1.0.jar
RUN mkdir -p /TMTIntegrator-v3.1.0/results && mv TMTIntegrator-v3.1.0.jar /TMTIntegrator-v3.1.0/

RUN curl -LOk https://github.com/Nesvilab/IonQuant/releases/download/1.7.17/IonQuant-1.7.17.jar
RUN mkdir /IonQuant-1.7.17 && mv IonQuant-1.7.17.jar /IonQuant-1.7.17/
RUN cp -r /MSFragger-3.4/ext /IonQuant-1.7.17/

COPY src/collect_psm.sh /usr/local/bin/
COPY src/package.sh /usr/local/bin/
COPY src/mv.sh /etc/profile.d/

CMD ["philosopher"]