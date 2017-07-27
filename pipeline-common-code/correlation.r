
## install and load libraries automatically
# library (psych)
if (!require("pacman")) install.packages ("pacman")
pacman::p_load (psych)


args <- commandArgs (TRUE)
pvalue <- 0.05


# id 1-jid.max are for calculating correlations
# id 0 is for assembing results and plotting
jid <- as.numeric (args[1])
jid.max <- as.numeric (args[2])
if (jid > jid.max) stop ("Invalid job id")

prefix <- toString (args[3])
pome.file <- toString (args[4])


read.matrix <- function (f) {
  d <- read.csv (f)
  rownames (d) <- d[,'GeneSymbol']
  return (d)
}



if (jid > 0) {
  pome <- read.matrix (pome.file)
  pome.data <- t (pome[,-1])
  
  n <- ncol (pome.data)
  job.size <- ceiling (n / jid.max)
  index.start <- (jid-1) * job.size + 1
  index.end <- min (jid * job.size, n)
  
  
  data.subset <- pome.data [, index.start:index.end]
  pome.corr <- corr.test (pome.data, data.subset, method='pearson', adjust='none')
  
  write.csv (round (pome.corr$r, digits=4), paste (prefix, '-output/pome-corr', jid, '.csv', sep=''))
  write.csv (round (pome.corr$p, digits=4), paste (prefix, '-output/pome-pval', jid, '.csv', sep=''))
}



if (jid == 0) {
  # read tables
  cat ('Reading tables ...\n')
  pome.corr <- NULL
  for (i in 1:jid.max) {
    if (i %% 10 == 0) cat (' job id:', i, '\n')
    if (i==1) {
      pome.corr$r <- read.csv (paste (prefix, '-output/pome-corr', i, '.csv', sep=''), row.names=1)
      pome.corr$p <- read.csv (paste (prefix, '-output/pome-pval', i, '.csv', sep=''), row.names=1)
    } else {
      pome.corr$r <- cbind (pome.corr$r, read.csv (paste (prefix, '-output/pome-corr', i, '.csv', sep=''), row.names=1))
      pome.corr$p <- cbind (pome.corr$p, read.csv (paste (prefix, '-output/pome-pval', i, '.csv', sep=''), row.names=1))
    }
  }
  
  # adjust p value
  cat ('Adjusting p-values ...\n')
  for (i in 1:ncol(pome.corr$p)) {
    pome.corr$p[,i] <- p.adjust (pome.corr$p[,i], method='BH')
  }
  
  # write out tables
  cat ('Writing consolidated tables ...\n')
  write.csv (pome.corr$r, paste (prefix, '-corr.csv', sep=''))
  write.csv (pome.corr$p, paste (prefix, '-pval.csv', sep=''))
  
  
  # count number of significant event for each (CNA) gene
  sig_pome <- sign (pome.corr$r) * as.numeric (pome.corr$p < pvalue)
  sig_pome_cnt <- apply (sig_pome, 2, function (x) sum (x != 0))
    
  # save significant event counts
  geneloc <- read.csv ('gene-location.csv')
  pome_events <- merge (geneloc, data.frame (GeneID=colnames (sig_pome), SignificantEvents=sig_pome_cnt),
                        by.x='HGNCsymbol', by.y='GeneID')
  write.csv (pome_events, paste (prefix, '-sigevents.csv', sep=''), row.names=FALSE)
}


