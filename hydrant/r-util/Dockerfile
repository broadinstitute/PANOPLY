
FROM r-base
MAINTAINER manidr@broadinstitute.org

# external dependencies for installing pacman library in R
# with pacman, any missing libraries will be automatically added
RUN apt-get update
RUN apt-get -t unstable install -y libssl-dev
RUN apt-get -t unstable install -y libcurl4-openssl-dev

# install packages
RUN echo 'install.packages(c( "Cairo", "ComplexHeatmap",     "MASS",     "MethComp",     "NMF",     "PerformanceAnalytics",     "RColorBrewer",     "RankAggreg",     "RobustRankAggreg",     "ape",     "bpca",     "caret",     "circlize",    "cluster",     "dplyr",    "e1071",     "fastcluster",     "ggplot2",     "glmnet",     "gplots",     "lattice",     "limma",     "lme4",     "maptools",     "mclust",     "misc3d",     "mixtools",     "nlme",     "optparse",    "pacman",    "pamr",     "parmigene",     "psych",     "randomForest",     "reshape",     "rgl",     "rhdf5",     "samr",     "scales",     "scatterplot3d",     "smacof",     "sn",     "tensor",     "tools",     "verification" ), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R  \
   && Rscript /tmp/packages.R

# bioconductor libraries
RUN echo 'source("https://bioconductor.org/biocLite.R"); biocLite("Biobase"); biocLite("graph"); biocLite("Rgraphviz"); biocLite("impute"); biocLite("rhdf5")' > /tmp/biocpkgs.R \
   && Rscript /tmp/biocpkgs.R

# bioconductor libraries for WGCNA (start by updating BiocGenerics)
RUN echo 'source("https://bioconductor.org/biocLite.R"); biocLite("BiocGenerics"); biocLite("preprocessCore"); biocLite("GO.db"); biocLite("AnnotationDbi")' > /tmp/biocpkgs-wgcna.R \
   && Rscript /tmp/biocpkgs-wgcna.R
# WGCNA for fast correlations
RUN echo 'install.packages("WGCNA")' > /tmp/packages-wgcna.R \
   && Rscript /tmp/packages-wgcna.R

# estimate library
RUN echo 'install.packages ("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)' > /tmp/packages-estimate.R \
   && Rscript /tmp/packages-estimate.R

COPY R-utilities/ /prot/proteomics/Projects/R-utilities

