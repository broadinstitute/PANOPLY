
## Configuration file for CMAP analysis
## Values set here supercede the defaults 
##  (defaults are used if this config file is not provided as an input)


# CMAP setup/initialization
cna.threshold <- 0.3               # copy number up/down if abs (log2(copy number) - 1) is > cna.threshold (default: 0.3)
cna.effects.threshold <- 15        # min number of samples with up/down copy number to include gene for CMAP analysis (default: 15)
min.sigevents <- 20                # gene must have at least this many significant trans events to be considered (default: 20)  
max.sigevents <- 1800              # if a gene has more then max.sigevents, the top max.sigevents will be chosen
                                   #  (having >1999 items in a geneset results in a crash -- unix buffer limits?)
top.N <- 500                       # maximum number of genes to run CMAP enrichment on (default: 500)
fdr.pvalue <- 0.05                 # FDR for CNA correlations (default: 0.05)
log.transform <-FALSE		   # if TRUE, log transform input data
must.include.genes <- NULL         # genes that must be included in the CMAP analysis (vector, default: NULL)

# CMAP connectivity score and statistical significance
cis.fdr <- fdr.pvalue              # FDR for cis-correlations (default: 0.05)
legacy.score <- FALSE              # if TRUE, legacy connectivity score will be calculated (using mean rank points), with permutation FDR
                                   # if FALSE, enrichement will be based on fisher test of outlier scores, with BH-FDR (default)
rankpt.n <- 4                      # number of CMAP profiles to consider when calculating mean rank point (default: 4)
mean.rankpt.threshold <- 85        # min value of mean rank point for gene signature to be considered enriched (default: 85)
cmap.fdr <- 0.25                   # BH-FDR threshold for fisher test of outlier scores, for gene to be considered enriched

# CMAP annotation
alpha <- fdr.pvalue		   # p-value threshold for cmap profile zscores and enrichments
