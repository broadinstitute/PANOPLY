FROM r-base:3.6.3
MAINTAINER proteogenomics@broadinstitute.org

RUN apt-get update
RUN apt-get -t unstable install -y libssl-dev
RUN apt-get -t unstable install -y libcurl4-openssl-dev
RUN apt-get install -y jags

COPY packages /prot/proteomics/Projects/PGDAC/packages

RUN echo 'if( !requireNamespace( "BiocManager", quietly = TRUE ) ) \
  install.packages( "BiocManager", \
    repos = "http://cran.us.r-project.org", dependencies = TRUE ); \
  BiocManager::install( version = "3.10" ); \
  BiocManager::install( c( "Biobase", \
	"graph", \
	"Rgraphviz", \
	"impute", \
	"BiocGenerics", \
	"preprocessCore", \
	"GO.db", \
        "rhdf5", \
        "prada", \
	"AnnotationDbi", \
        "ComplexHeatmap", \
        "SummarizedExperiment" ) );' > /tmp/biocpkgs.R \
  && Rscript /tmp/biocpkgs.R

RUN echo 'install.packages ( "estimate", \
	repos = "http://r-forge.r-project.org", dependencies = TRUE )' \
	> /tmp/packages-estimate.R \
  && Rscript /tmp/packages-estimate.R

RUN echo 'install.packages( "pacman", repos = "http://cran.us.r-project.org", \
        dependencies = TRUE ); library( pacman ); \
	p_install( data.table ); \
	p_install( MASS ); \
	p_install( NMF ); \
	p_install( PerformanceAnalytics ); \
	p_install( RColorBrewer ); \
	p_install( RankAggreg ); \
	p_install( RobustRankAggreg ); \
	p_install( ape ); \
	p_install( caret ); \
	p_install( circlize ); \
	p_install( cluster ); \
	p_install( e1071 ); \
	p_install( fastcluster ); \
	p_install( glmnet ); \
	p_install( gplots ); \
	p_install( lattice ); \
	p_install( lme4 ); \
	p_install( maptools ); \
	p_install( mclust ); \
	p_install( misc3d ); \
	p_install( mixtools ); \
	p_install( nlme ); \
	p_install( optparse ); \
	p_install( parmigene ); \
	p_install( psych ); \
	p_install( randomForest ); \
	p_install( reshape ); \
	p_install( samr ); \
	p_install( scales ); \
	p_install( scatterplot3d ); \
	p_install( smacof ); \
	p_install( sn ); \
	p_install( tensor ); \
	p_install( tidyverse ); \
	p_install( tools ); \
	p_install( verification ); \
	p_install( WGCNA ); \
	p_install( factoextra ); \
	p_install( glue ); \
	p_install( magrittr ); \
	p_install( doParallel ); \
	p_install( foreach ); \
	p_install( pheatmap ); \
	p_install( NbClust ); \
	p_install( rjson ); \
	p_install( GlobalOptions); \
	p_install( GetoptLong ); \
        p_install( coda ); \
        p_install( rjags ); \
        p_install( limma );' > /tmp/packages.R  && Rscript /tmp/packages.R


RUN echo 'install.packages( \
	"/prot/proteomics/Projects/PGDAC/packages/cmapR_0.99.19.tar.gz", \
	type="source", repos=NULL );' > /tmp/packages.R && Rscript /tmp/packages.R

RUN echo 'install.packages( \
        "/prot/proteomics/Projects/PGDAC/packages/MethComp_1.30.0.tar.gz", \
        type="source", repos=NULL );' > /tmp/packages.R && Rscript /tmp/packages.R

RUN apt-get -t unstable install -y gfortran -std=legacy
RUN echo 'install.packages( \
	"/prot/proteomics/Projects/PGDAC/packages/pamr_1.56.1.tar.gz", \
	type="source", repos=NULL );' > /tmp/packages.R && Rscript /tmp/packages.R

## Check which packages got installed and which failed
COPY src/test_packages.r /tmp/test_packages.r
RUN  Rscript /tmp/test_packages.r > /tmp/check_pkgs_log.txt
