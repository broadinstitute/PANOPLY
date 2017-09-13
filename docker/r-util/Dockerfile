
FROM r-base
MAINTAINER manidr@broadinstitute.org

# external dependencies for installing pacman library in R
# with pacman, any missing libraries will be automatically added
RUN apt-get update
RUN apt-get -t unstable install -y libssl-dev
RUN apt-get -t unstable install -y libcurl4-openssl-dev

# install packages
RUN echo 'install.packages(c(     "MASS",     "MethComp",     "NMF",     "PerformanceAnalytics",     "RColorBrewer",     "RankAggreg",     "RobustRankAggreg",     "bpca",     "caret",     "e1071",     "fastcluster",     "ggplot2",     "glmnet",     "gplots",     "lattice",     "limma",     "lme4",     "maptools",     "mclust",     "misc3d",     "mixtools",     "nlme",     "pacman",    "pamr",     "parmigene",     "psych",     "randomForest",     "reshape",     "rgl",     "rhdf5",     "samr",     "scales",     "scatterplot3d",     "smacof",     "sn",     "tensor",     "tools",     "verification" ), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R  \
   && Rscript /tmp/packages.R
# bioconductor libraries
RUN echo 'source("https://bioconductor.org/biocLite.R"); biocLite("Biobase"); biocLite("graph"); biocLite("Rgraphviz"); biocLite("impute"); biocLite("rhdf5")' > /tmp/biocpkgs.R \
   && Rscript /tmp/biocpkgs.R

COPY R-utilities/ /prot/proteomics/Projects/R-utilities

