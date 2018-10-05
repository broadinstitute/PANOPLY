


## generate input data files for CMAP analysis
generate.cmap.input <- function (target.cna.dir,
                                 cmap.kd.genelist, 
                                 input.genes='cmap-input-genes-list.txt',   # genes listed here will always be included
                                 group='all',                # subgroup name on which CNA analysis was run
                                 dtype='pome',               # dtype is either pome or mrna
                                 generate.permuted.genesets=0,  # number of permuted datasets to generate
                                                             # if >0, generate permuted genesets; else output actual genesets
                                 cna.threshold=0.3,          # copy number up/down if abs (log2(copy number) - 1) is > 0.3
                                 cna.effects.threshold=15,   # min number of tumors with up/down copy number to include gene for CMAP analysis
                                 min.sigevents=20,           # gene must have at least this many significant trans events to be considered  
                                 max.sigevents=1800,         # if a gene has more then max.sigevents, the top max.sigevents will be chosen
                                                             #  (having >1999 items in a geneset results in a crash -- unix buffer limits?)
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
  sigevents.all <- read.csv (file.path (target.cna.dir, paste (group, dtype, 'vs-cna-sigevents.csv', sep='-')))
  pvals <- read.csv (file.path (target.cna.dir, paste (group, dtype, 'vs-cna-pval.csv', sep='-')), row.names=1)
  data <- read.csv (file.path (target.cna.dir, paste (group, dtype, 'matrix.csv', sep='-')))

  # calculate cis correlation to cna -- this is needed to filter enriched signatures in connectivity.r
  cis.corr <- calc.cis.correlation (cna, data)
  all.genes <- as.character (cis.corr [,'gene'])
  write.csv (cis.corr, paste (group, dtype, 'cis-correlation.csv', sep='-'), row.names=FALSE, quote=FALSE)
  
  # filter to determine genes to run CMAP on
  sigevents <- sigevents.all
  keep <- sigevents[,'HGNCsymbol'] %in% selected.cna.genes & sigevents[,'SignificantEvents'] >= min.sigevents
  sigevents <- sigevents [keep,]
  event.count <- order (sigevents [,'SignificantEvents'], decreasing=TRUE)
  sigevents <- sigevents [ event.count, ]
  # top genes are picked based on (i) number of significant events; and
  # (ii) whether they have been knocked down in the CMAP profiles
  # if provided, genes listed in input.genes are always retained
  if (file.exists(input.genes)) {
    must.keep <- as.character (read.delim (input.genes, header=FALSE)[,1])
    sigevents <- rbind (sigevents.all [sigevents.all [,'HGNCsymbol'] %in% must.keep, ],
                        sigevents)
  }
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
    
  # output files for genesets
  cmap.updn.genes.file <- paste (group, 'cmap', dtype, 'updn-genes.gmt', sep='-')
  permuted.genes.prefix <- paste (group, 'cmap', dtype, 'permuted-genes', sep='-')
  if (file.exists (cmap.updn.genes.file)) file.remove (cmap.updn.genes.file)  # genesets are appended to file -- start fresh
  if (generate.permuted.genesets > 0) {
    for (p in 1:generate.permuted.genesets) 
      if (file.exists (sprintf ("%s-%03d.gmt", permuted.genes.prefix, p))) 
        file.remove (sprintf ("%s-%03d.gmt", permuted.genes.prefix, p))
  }
  
  write.geneset <- function (gene, up.list, dn.list, output.file) {
    # function to write out a single geneset to a specified file
    cat (gene, '_CNA', v, '\t', 'na\t', paste (up.list, collapse=';u\t'), ';u\t', paste (dn.list, collapse=';d\t'), ';d', '\n',
         file=output.file, append=TRUE, sep='')
  }
  
  for (g in top.genes) {
    cat ('  gene', g, '\n')
    
    # get CNA "class vector"
    g.cna <- cna.data [ which (cna.genes == g), ]
    
    # get trans genes
    g.pvals <- pvals [, g]
    max.sigevent.pval <- g.pvals [ order (g.pvals) ] [max.sigevents]
    siggenes <- rownames (pvals) [g.pvals < min (fdr.pvalue, max.sigevent.pval, na.rm=TRUE)]
    
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
          write.geneset (g, up, dn, cmap.updn.genes.file)
          if (generate.permuted.genesets > 0) {
            ## generate permuted datasets for FDR calculation
            for (p in 1:generate.permuted.genesets) {
              # for each gene, generate a random set of trans genes
              # (generate up and dn trans genes together to prevent overlap)
              trans.p <- sample (all.genes, length (c (up, dn)))
              up.p <- trans.p [ 1:length(up) ]
              dn.p <- rev (trans.p) [ 1:length(dn) ]
              write.geneset (g, up.p, dn.p, sprintf ("%s-%03d.gmt", permuted.genes.prefix, p))
            }
          }
        }
      }
    }
  }
}




calc.cis.correlation <- function (data1, data2) {
  ## support function to calculate cis-correlation
  # data prep
  rownames (data1) <- data1 [,'GeneSymbol']
  rownames (data2) <- data2 [,'GeneSymbol']
  common.genes <- intersect (data1[,'GeneSymbol'], data2[,'GeneSymbol'])
  data1 <- t (data1 [common.genes,-1])
  data2 <- t (data2 [common.genes,-1])
  
  corr <- sapply (1:length(common.genes),
                  function (x) {
                    result <- cor.test (data1[,x], data2[,x], method='pearson')
                    return (c (result$estimate, result$p.value))
                  })
  retval <- cbind (common.genes, t (corr), p.adjust (corr[2,], method='BH'))
  colnames (retval) <- c ('gene', 'correlation', 'pvalue', 'adj.pvalue')
  return (retval)
}




## read argument options from command line
args <- commandArgs (TRUE)

cna.dir <- toString (args[1])
kd.list <- toString (args[2])
gr <- toString (args[3])
typ <- toString (args[4])
perm <- as.numeric (args[5])
logtr <- ifelse (toString (args[6])=='TRUE', TRUE, FALSE)

## execute generator function
generate.cmap.input (target.cna.dir=cna.dir, cmap.kd.genelist=kd.list, group=gr,
                     dtype=typ, generate.permuted.genesets=perm, log.transform=logtr)

