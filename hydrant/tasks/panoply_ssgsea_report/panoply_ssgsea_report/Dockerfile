FROM broadcptacdev/panoply_rmd_report:latest

LABEL maintainer="proteogenomics@broadinstitute.org"
LABEL description="Docker for generating Rmarkdown reports for different PGDAC modules."

ENV wd /home/pgdac/

COPY . $wd
WORKDIR $wd


RUN R -e "pacman::p_install('pheatmap')"
RUN R -e "pacman::p_install('ggpubr')"
RUN R -e "pacman::p_install('ggrepel')"
RUN R -e "pacman::p_install('RColorBrewer')"
RUN R -e "pacman::p_install('caTools')"
RUN R -e "pacman::p_install('TSP')"
RUN R -e "pacman::p_install('seriation')"





