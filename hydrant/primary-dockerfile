FROM broadcptac/r_util:2
MAINTAINER rkothadi@broadinstitute.org

COPY src /prot/proteomics/Projects/PGDAC/src
COPY data /prot/proteomics/Projects/PGDAC/data
COPY packages /prot/proteomics/Projects/PGDAC/packages

RUN R -e "library( pacman ); p_install( 'factoextra'     )"
RUN R -e "library( pacman ); p_install( 'glue'           )"
RUN R -e "library( pacman ); p_install( 'magrittr'       )"
RUN R -e "library( pacman ); p_install( 'doParallel'     )"
RUN R -e "library( pacman ); p_install( 'foreach'        )"
RUN R -e "library( pacman ); p_install( 'pheatmap'       )"
RUN R -e "library( pacman ); p_install( 'RcolorBrewer'   )"
RUN R -e "library( pacman ); p_install( 'NbClust'        )"
RUN R -e "library( pacman ); p_install( 'rjson'          )"
RUN R -e "library( pacman ); p_install( 'GlobalOptions'  )"
RUN R -e "library( pacman ); p_install( 'GetoptLong'     )"
RUN R -e "library( pacman ); p_install( 'ComplexHeatmap' )"
RUN R -e "library( pacman ); p_install( 'dplyr'          )"
