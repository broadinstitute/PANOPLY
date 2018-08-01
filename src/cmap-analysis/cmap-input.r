


## generate input data files for CMAP analysis
generate.cmap.input <- function (target.cna.dir,
                                 cmap.kd.genelist, 
                                 group='all',                # subgroup name on which CNA analysis was run
                                 dtype='pome',               # dtype is either pome or mrna
                                 generate.permuted.genesets=0,  ### NEEDS IMPLEMENTATION
                                                             # if >0, generate permuted genesets; else output actual genesets
                                 cna.threshold=0.3,          # copy number up/down if abs (log2(copy number) - 1) is > 0.3
                                 cna.effects.threshold=15,   # min number of tumors with up/down copy number to include gene for CMAP analysis
                                 min.sigevents=20,           # gene must have at least this many significant trans events to be considered  
                                 top.N=500,
                                 fdr.pvalue=0.05,
                                 exclude.cis=FALSE,
                                 log.transform=FALSE)
{
  ## read CNA data and log transform
  cna <- read.csv (file.path (target.cna.dir, sprintf ('%s-cna-matrix.csv', group)))
  cna.genes <- as.character (cna [,1])

  cna.data <- cna [,2:ncol(cna)]
  if (log.transform) cna.data <- log2 (cna.data) - 1
  # convert to up (+1) / down (-1) calls
  cna.data <- sign (cna.data) * (abs (cna.data) > cna.threshold)

  updn.count <- apply (cna.data, 1, function (x) sum (x != 0))
  cna.keep <- updn.count >= cna.effects.threshold
  selected.cna.genes <- cna.genes [cna.keep]


  # read in all required tables
  sigevents <- read.csv (file.path (target.cna.dir, paste (group, dtype, 'vs-cna-sigevents.csv', sep='-')))
  pvals <- read.csv (file.path (target.cna.dir, paste (group, dtype, 'vs-cna-pval.csv', sep='-')), row.names=1)
  data <- read.csv (file.path (target.cna.dir, paste (group, dtype, 'matrix.csv', sep='-')))

  # filter to determine genes to run CMAP on
  keep <- sigevents[,'HGNCsymbol'] %in% selected.cna.genes & sigevents[,'SignificantEvents'] >= min.sigevents
  sigevents <- sigevents [keep,]
  event.count <- order (sigevents [,'SignificantEvents'], decreasing=TRUE)
  sigevents <- sigevents [ event.count, ]
  # top genes are picked based on (i) number of significant events; and
  # (ii) whether they have been knocked down in the CMAP profiles
  cmap.sh.genes <- as.character (read.delim (cmap.kd.genelist, header=FALSE)[,1])
  top.genes <- NULL
  i <- 1
  while (length(top.genes) <= top.N && i <= nrow(sigevents)) {
    g <- as.character (sigevents [i, 'HGNCsymbol'])
    if (g %in% cmap.sh.genes) top.genes <- c (top.genes, g)
    i <- i+1
  }
  write.csv (sigevents [ sigevents[,'HGNCsymbol'] %in% top.genes, ],
             paste (group, 'cmap', dtype, 'gene-list.csv', sep='-'), row.names=FALSE)
    
  cmap.updn.genes.file <- paste (group, 'cmap', dtype, 'updn-genes.gmt', sep='-')
  for (g in top.genes) {
    cat ('  gene', g, '\n')
    
    # get CNA "class vector"
    g.cna <- cna.data [ which (cna.genes == g), ]
    
    # get trans genes
    g.pvals <- pvals [, g]
    siggenes <- rownames (pvals) [g.pvals < fdr.pvalue]
    
    # get expression data for trans genes
    g.alldata <- data [ data[,'GeneSymbol'] %in% siggenes,]
    g.data <- g.alldata [,-1]
    g.genes <- g.alldata [,1]
    
    # calculate median expression for CNA del (-1), neutral (0) and amp (+1)
    g.median <- t ( apply (g.data, 1,
                           function (x) c (del=median (x[g.cna==-1], na.rm=TRUE),
                                           neu=median (x[g.cna==0], na.rm=TRUE),
                                           amp=median (x[g.cna==1], na.rm=TRUE))) )
    
    # for CNA-DEL and CNA-AMP, determine UP/DN genes (by comparing to CNA-NEU)
    g.updn <- t ( apply (g.median, 1, function (x) c (x['del'] > x['neu'],
                                                      x['amp'] > x['neu'])) )
    
    # write CMAP up / dn gmt files
    # (using gene names -- convert to affy ids later -- had to install mygene package)
    for (v in c ('del', 'amp')) {
      if (v %in% colnames (g.updn)) {
        # there may be no amp and/or del samples
        up <- g.genes [g.updn[,v]]
        dn <- g.genes [!g.updn[,v]]
        # also exclude any missing values, which can happen when del/neu/amp is NA
        up <- up [ !is.na (up) ]
        dn <- dn [ !is.na (dn) ]
        if (exclude.cis) {
          # exclude the (cis) gene itself
          up <- setdiff (up, g)
          dn <- setdiff (dn, g)
        }
        if (length(up) > 0 && length(dn) > 0) {
          cat (g, '_CNA', v, '\t', 'na\t', paste (up, collapse=';u\t'), ';u\t', paste (dn, collapse=';d\t'), ';d', '\n',
               file=cmap.updn.genes.file, append=TRUE, sep='')
        }
      }
    }
  }
}


## read argument options from command line
args <- commandArgs (TRUE)

cna.dir <- toString (args[1])
kd.list <- toString (args[2])
gr <- toString (args[3])
typ <- toString (args[4])
perm <- as.numeric (args[5])

## execute generator function
generate.cmap.input (target.cna.dir=cna.dir, cmap.kd.genelist=kd.list, group=gr,
                     dtype=typ, generate.permuted.genesets=perm)

