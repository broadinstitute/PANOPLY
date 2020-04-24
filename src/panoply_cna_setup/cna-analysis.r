
## install and load libraries automatically
# library (psych)
if (!require("pacman")) install.packages ("pacman")
pacman::p_load (WGCNA)
allowWGCNAThreads()

# NB: WGCNA is a fast correlation calculation library -- 2-3 orders of magnitude faster than previous library (psych).
# Automatic installation may fail; the following may need to be
# installed manually: BioGenerics, preprocessCore, GO.db, AnnotationDbi

source ('generate-cna-plots.r')
source ('config.r')

args <- commandArgs (TRUE)
#pvalue <- 0.05
pvalue <- fdr_cna_corr

# id 1-jid.max are for calculating correlations
# id 0 is for assembing results and plotting
jid <- as.numeric (args[1])
jid.max <- as.numeric (args[2])
if (jid > jid.max) stop ("Invalid job id")

prefix <- toString (args[3])
mrna.file <- toString (args[4])
cna.file <- toString (args[5])
pome.file <- toString (args[6])


read.matrix <- function (f) {
  d <- read.csv (f)
  rownames (d) <- d[,'GeneSymbol']
  return (d)
}


if (jid > 0) {
  mrna <- read.matrix (mrna.file)
  cna <- read.matrix (cna.file)
  pome <- read.matrix (pome.file)
  
  ### Do NOT subset to common genes here -- that restricts correlations to common genes only
  ### instead, calculate all correlations and subset to common genes before plotting
  ##
  ## common.genes <- intersect (mrna[,'GeneSymbol'],
  ##                            intersect (pome[,'GeneSymbol'], cna[,'GeneSymbol']))
  ## mrna <- t (mrna [common.genes, -1])
  ## cna <- t (cna [common.genes, -1])
  ## pome <- t (pome [common.genes, -1])
  
  
  ## job.size <- ceiling (length (common.genes) / jid.max)
  ## index.start <- (jid-1) * job.size + 1
  ## index.end <- min (jid * job.size, length (common.genes))

  mrna <- t (mrna[,-1])
  cna <- t (cna[,-1])
  pome <- t (pome[,-1])

  # parallellize by splitting the CNA table into strips
  job.size <- ceiling (ncol (cna) / jid.max)
  index.start <- (jid-1) * job.size + 1
  index.end <- min (jid * job.size, ncol(cna))



  cna.subset <- cna [, index.start:index.end]
  mrna.vs.cna <- corAndPvalue (mrna, cna.subset)
  pome.vs.cna <- corAndPvalue (pome, cna.subset)
  
  write.csv (round (mrna.vs.cna$cor, digits=4), paste (prefix, '-output/mrna-vs-cna-corr', jid, '.csv', sep=''))
  write.csv (round (mrna.vs.cna$p, digits=4), paste (prefix, '-output/mrna-vs-cna-pval', jid, '.csv', sep=''))
  write.csv (round (pome.vs.cna$cor, digits=4), paste (prefix, '-output/pome-vs-cna-corr', jid, '.csv', sep=''))
  write.csv (round (pome.vs.cna$p, digits=4), paste (prefix, '-output/pome-vs-cna-pval', jid, '.csv', sep=''))
}


if (jid == 0) {
  # read tables
  cat ('Reading tables ...\n')
  mrna.vs.cna <- pome.vs.cna <- NULL

  # check if consolidated tables are present -- many times, the correlation calcs run, but plot
  # generation fails (esp when using gold); in such cases, read the consolidated tables and proceed
  mrna.corr.file <- paste (prefix, '-mrna-vs-cna-corr.csv', sep='')
  mrna.pval.file <- paste (prefix, '-mrna-vs-cna-pval.csv', sep='')
  pome.corr.file <- paste (prefix, '-pome-vs-cna-corr.csv', sep='')
  pome.pval.file <- paste (prefix, '-pome-vs-cna-pval.csv', sep='')

  if ( all (file.exists (c (mrna.corr.file, mrna.pval.file, pome.corr.file, pome.pval.file))) ) {
    # files exist -- just read and proceed
    cat ('Conslidated tables exist ... \n')
    cat ('Reading consolidated tables ...\n')
    mrna.vs.cna$r <- read.csv (mrna.corr.file, row.names=1)
    mrna.vs.cna$p <- read.csv (mrna.pval.file, row.names=1)
    pome.vs.cna$r <- read.csv (pome.corr.file, row.names=1)
    pome.vs.cna$p <- read.csv (pome.pval.file, row.names=1)
  } else {
    # files do not exist -- create them by assembling results from parallel runs
    for (i in 1:jid.max) {
      if (i %% 10 == 0) cat (' job id:', i, '\n')
      if (i==1) {
        mrna.vs.cna$r <- read.csv (paste (prefix, '-output/mrna-vs-cna-corr', i, '.csv', sep=''), row.names=1)
        mrna.vs.cna$p <- read.csv (paste (prefix, '-output/mrna-vs-cna-pval', i, '.csv', sep=''), row.names=1)
        pome.vs.cna$r <- read.csv (paste (prefix, '-output/pome-vs-cna-corr', i, '.csv', sep=''), row.names=1)
        pome.vs.cna$p <- read.csv (paste (prefix, '-output/pome-vs-cna-pval', i, '.csv', sep=''), row.names=1)
      } else {
        mrna.vs.cna$r <- cbind (mrna.vs.cna$r, read.csv (paste (prefix, '-output/mrna-vs-cna-corr', i, '.csv', sep=''), row.names=1))
        mrna.vs.cna$p <- cbind (mrna.vs.cna$p, read.csv (paste (prefix, '-output/mrna-vs-cna-pval', i, '.csv', sep=''), row.names=1))
        pome.vs.cna$r <- cbind (pome.vs.cna$r, read.csv (paste (prefix, '-output/pome-vs-cna-corr', i, '.csv', sep=''), row.names=1))
        pome.vs.cna$p <- cbind (pome.vs.cna$p, read.csv (paste (prefix, '-output/pome-vs-cna-pval', i, '.csv', sep=''), row.names=1))
      
      }
    }
    

    ## subset to common genes so that the CNA plot will have the same x and y axes
    ## full results (without subsetting) is contained in the output/table
    common.genes <- intersect (intersect (rownames (mrna.vs.cna$r), colnames (mrna.vs.cna$r)),
                               intersect (rownames (pome.vs.cna$r), colnames (pome.vs.cna$r)))
    mrna.vs.cna$r <- mrna.vs.cna$r [common.genes, common.genes]
    mrna.vs.cna$p <- mrna.vs.cna$p [common.genes, common.genes]
    pome.vs.cna$r <- pome.vs.cna$r [common.genes, common.genes]
    pome.vs.cna$p <- pome.vs.cna$p [common.genes, common.genes]

  
    # adjust p value
    cat ('Adjusting p-values ...\n')
    for (i in 1:ncol(mrna.vs.cna$p)) mrna.vs.cna$p[,i] <- p.adjust (mrna.vs.cna$p[,i], method='BH')
    for (i in 1:ncol(pome.vs.cna$p)) pome.vs.cna$p[,i] <- p.adjust (pome.vs.cna$p[,i], method='BH')
    
    
    # write out tables
    cat ('Writing consolidated tables ...\n')
    write.csv (mrna.vs.cna$r, mrna.corr.file)
    write.csv (mrna.vs.cna$p, mrna.pval.file)
    write.csv (pome.vs.cna$r, pome.corr.file)
    write.csv (pome.vs.cna$p, pome.pval.file)
  }
  
  
  # create final tables
  cna_mrna <- as.matrix (sign (mrna.vs.cna$r) * as.numeric (mrna.vs.cna$p < pvalue))
  cna_pome <- as.matrix (sign (pome.vs.cna$r) * as.numeric (pome.vs.cna$p < pvalue))

  
  # count number of significant event for each (CNA) gene
  cna_mrna_cnt <- apply (cna_mrna, 2, function (x) sum (x != 0))
  cna_pome_cnt <- apply (cna_pome, 2, function (x) sum (x != 0))
  
  # read in tables needed for plotting
  geneloc <- read.csv ('gene-location.csv')
  chrlength <- read.csv ('chr-length.csv', header=FALSE)
  
  # save significant event counts
  cna_in_mrna <- match (colnames (cna_mrna), rownames (cna_mrna))
  mrna_events <- merge (geneloc, data.frame (GeneID=colnames (cna_mrna), 
                                             SignificantEvents=cna_mrna_cnt,
                                             SignificantCisEffect=(diag (cna_mrna [cna_in_mrna, ]) != 0)),
                        by.x='HGNCsymbol', by.y='GeneID')
  cna_in_pome <- match (colnames (cna_pome), rownames (cna_pome))
  pome_events <- merge (geneloc, data.frame (GeneID=colnames (cna_pome), 
                                             SignificantEvents=cna_pome_cnt,
                                             SignificantCisEffect=(diag (cna_pome [cna_in_pome, ]) != 0)),
                        by.x='HGNCsymbol', by.y='GeneID')
  write.csv (mrna_events, paste (prefix, '-mrna-vs-cna-sigevents.csv', sep=''), row.names=FALSE)
  write.csv (pome_events, paste (prefix, '-pome-vs-cna-sigevents.csv', sep=''), row.names=FALSE)
  
  
  # call plotting routine
  cat ('Creating plot ...\n')
  Plot_cis_trans_effect (cna_mrna, cna_pome, geneloc, chrlength, paste (prefix, '-cna-plot.png', sep=''))
}

