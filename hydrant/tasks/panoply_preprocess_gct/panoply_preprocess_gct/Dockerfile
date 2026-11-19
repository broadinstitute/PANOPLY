FROM broadcptacdev/panoply_ssgsea:latest
MAINTAINER proteogenomics@broadinstitute.org

## Copy SRC Files into Docker
ENV wd=/home/pgdac/
COPY . $wd
WORKDIR $wd

# Annotation packages used to map RefSeq accessions / gene symbols to UniProt ids
# Note: the four org.*.eg.db packages are large (~1.3 GB together); consider dropping the organisms you don't need.
RUN Rscript -e 'options(Ncpus = parallel::detectCores()); \
      install.packages(c("AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dr.eg.db"), repos = c(getOption("repos"), "https://p3m.dev/bioconductor/__linux__/noble/2026-06-23/packages/3.23/bioc", "https://p3m.dev/bioconductor/__linux__/noble/2026-06-23/packages/3.23/data/annotation")); \
      stopifnot(sapply(c("AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dr.eg.db"), requireNamespace, quietly = TRUE))' # stop if packages didn't build

## cleanup files from package installs
RUN rm -rf /tmp/downloaded_packages /tmp/Rtmp*
