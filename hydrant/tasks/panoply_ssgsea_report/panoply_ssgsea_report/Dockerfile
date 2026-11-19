FROM broadcptacdev/panoply_ssgsea:latest

LABEL maintainer="proteogenomics@broadinstitute.org"
LABEL description="Docker for generating the Rmarkdown report of the panoply_ssgsea module."

# Packages for rendering the Rmarkdown report of the panoply_ssgsea_report module
#   pandoc: required by rmarkdown to convert reports to HTML
#   kableExtra: formatted tables in the report (rmarkdown and knitr are installed with tidyverse)
RUN apt-get update \
	&& apt-get install -y --no-install-recommends pandoc \
	&& rm -rf /var/lib/apt/lists/* \
	&& install2.r --error --skipinstalled -n $(nproc) kableExtra

## cleanup files from package installs
RUN rm -rf /tmp/downloaded_packages /tmp/Rtmp*

## Copy SRC Files into Docker
ENV wd=/home/pgdac/
COPY . $wd
WORKDIR $wd