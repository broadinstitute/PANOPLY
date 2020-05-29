install.packages( "pacman" )
library( pacman )
p_load( glue )


libs <- c( "Biobase", "graph", "Rgraphviz", "impute", "BiocGenerics", "preprocessCore", "GO.db", "rhdf5", "prada", "AnnotationDbi", "estimate", "data.table", "MASS", "NMF", "PerformanceAnalytics", "RColorBrewer", "RankAggreg", "RobustRankAggreg", "ape", "caret", "circlize", "cluster", "dplyr", "e1071", "fastcluster", "ggplot2", "glmnet", "gplots", "lattice", "lme4", "maptools", "mclust", "misc3d", "mixtools", "nlme", "optparse", "parmigene", "psych", "randomForest", "reshape", "samr", "scales", "scatterplot3d", "smacof", "sn", "tensor", "tools", "verification", "WGCNA", "factoextra", "glue", "magrittr", "doParallel", "foreach", "pheatmap", "NbClust", "rjson", "GlobalOptions", "GetoptLong", "MethComp", "limma" )

libs <- c( libs, "ComplexHeatmap", "pamr", "cmapR" )

ipkg <- installed.packages()
fpkg <- setdiff( libs, ipkg )
spkg <- intersect( libs, ipkg )
	
nos <- length( spkg )
ver <- c()
for ( num in 1:nos ){
  ver <- paste0( packageVersion( spkg[num] ), collapse = '.' )
  print( glue( "{spkg[num]} : {ver}\n" ) )
}

print( glue( "\n\n=====================================" ) )
print( glue( "\nFailed to install: " ) )
print( glue( "{paste0( fpkg, collapse = ', ' )}" ) )
print( glue( "\nTried to install: {length( libs )}" ) )
print( glue( "\nTotal installed: {nos}" ) )
