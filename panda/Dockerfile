FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:2.0.5

# install AnnotationDbi for protein keytype ID
RUN R -e "BiocManager::install(c('cmapR'))"
RUN R -e "BiocManager::install(c('AnnotationDbi', 'org.Hs.eg.db', 'EnsDb.Hsapiens.v79'))"
RUN R -e "install.packages('khroma')"
COPY R-utilities /prot/proteomics/Projects/R-utilities
COPY notebook_startup.sh /tmp/notebook_startup.sh
RUN bash /tmp/notebook_startup.sh
