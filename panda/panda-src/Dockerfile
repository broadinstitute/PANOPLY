FROM broadcptacdev/panda_config_libs:latest

COPY build-config.r /panda/build-config.r
COPY bin /panda/bin
COPY defaults /panda/defaults

# install AnnotationDbi for protein keytype ID
RUN R -e "BiocManager::install('AnnotationDbi'); BiocManager::install('org.Hs.eg.db')"
