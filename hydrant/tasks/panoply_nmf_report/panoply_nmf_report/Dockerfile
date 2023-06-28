FROM broadcptacdev/panoply_rmd_report:latest

# Install any libraries necessary to run the task
# RUN set -ex \
#     && apt-get update \
#     && apt-get install -y --no-install-recommends \
#         libcairo2-dev \
#         libxt-dev \
#     && rm -rf /var/lib/apt/lists/* \
#     && install2.r -e \
#         Cairo \
#         data.table \
#     && rm -rf /tmp/*
RUN apt-get update --allow-releaseinfo-change && apt-get install -y dos2unix

ENV wd /home/pgdac/
COPY . $wd
WORKDIR $wd
RUN dos2unix src/*