FROM broadcptacdev/panoply_common:latest
MAINTAINER proteogenomics@broadinstitute.org


## install webshot
RUN R -e "install.packages('webshot')"

## install and configure phantomJS (dependency of webshot)
RUN wget https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2
RUN tar -xf phantomjs-2.1.1-linux-x86_64.tar.bz2 -C /usr/local/share/  # untar package
RUN ln -sf /usr/local/share/phantomjs-2.1.1-linux-x86_64/bin/phantomjs /usr/local/bin  # create symbolic link to run phantomjs
RUN touch /tmp/openssl.cnf # create necessary config file
ENV OPENSSL_CONF /tmp/openssl.cnf # set relevant environmental variable for config file

COPY src /prot/proteomics/Projects/PGDAC/src


