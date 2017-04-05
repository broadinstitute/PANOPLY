


library (ggplot2)
library (reshape)

source ('preamble.r')
Source ('misc.r')
Source ('map-to-genes.r')



explore.correlations <- function (prefix, pome.gct.file, mrna.gct.file, pam50.cls, 
                                  exclude=NULL, FDR=0.01, correlation.type='pearson') {
  # explore gene (mrna) - protein correlation and sample-wise correlation between 
  # mrna expression and protein expression
  # any samples in the exclude vector (ala cls) are removed before correlation calcs/plots
  # n.b: "pome" is used to repesent proteome or phosphoproteome
  
  gi.gene.map <- read.delim (file.path (data.dir, protein.gene.map), sep='\t')

  # read cls and account for non-BC samples
  cls <- read.cls (pam50.cls)
  keep <- cls == "Basal" | cls == "Her2" | cls == "LumA" | cls == "LumB"
  cls <- cls [ keep ]
  if (!is.null (exclude)) exclude <- exclude [ keep ]

  
  read.gct.filter.cols <- function (gct.file, col1.name) {
    gdata <- read.gct (gct.file, check.names=FALSE)
    gdata <- gdata [, c (TRUE, TRUE, keep)]
    
    if (!is.null (exclude)) gdata <- gdata [, c (TRUE, TRUE, !exclude)]
    gdata <- gdata [, -1]
    colnames (gdata)[1] <- col1.name
    colnames (gdata)[-1] <- unlist (lapply (colnames (gdata)[-1],
                                            function (x) { strsplit (x, split='\\.')[[1]][1] }))
    colnames (gdata)[-1] <- sub ('\\.', '-', make.names (colnames (gdata)[-1], unique=TRUE) )
    return (gdata)
  }

  
  pome <- read.gct.filter.cols (pome.gct.file, protein.id.col)
  mrna <- read.gct.filter.cols (mrna.gct.file, gene.id.col)
  
  
  # reorder by PAM-50 class (cls file is identical for pome and mrna)
  # order (alphabetical): Basal, Her2, LumA, LumB
  # eliminate Normal samples since these don't have mRNA profiles
  if (!is.null (exclude)) cls <- cls [!exclude]
  
  new.order <- c (1, order (cls) + 1)
  pome <- pome [, new.order]
  mrna <- mrna [, new.order] 
  
  
  data <- merge (gi.gene.map, pome, by=protein.id.col)
  data <- merge (data, mrna, by=gene.id.col, suffixes=c ('.pome', '.mrna'))
  data <- process.duplicate.genes (data, genesym.col=1, data.cols=7:ncol(data), map.genes=FALSE)
  
  # data extends from start.col:end.col, with proteomics samples first, 
  # followed by mrna samples
  pome.cols <- grep ('\\.pome', colnames (data))
  mrna.cols <- grep ('\\.mrna', colnames (data))
  if ( length(pome.cols) != length(mrna.cols) ) stop (paste (type, 'and mRNA have different numbers of sample'))
  
  
  ##
  ## gene/protein correlation
  ##
  
  correlations <- apply (data, 1,
                         function (x) {
                           tryCatch (cor.test (as.numeric (x[pome.cols]), as.numeric (x[mrna.cols]), 
                                               method=correlation.type)$estimate,
                                     error = function (e) NA)
                         })
  pvalues <- apply (data, 1,
                    function (x) {
                      tryCatch (cor.test (as.numeric (x[pome.cols]), as.numeric (x[mrna.cols]), 
                                          method=correlation.type)$p.value,
                                error = function (e) NA)
                    })
  pvalues.fdr <- p.adjust (pvalues, method='fdr')
  
  data <- cbind (data, correlation=correlations, p.value=pvalues, p.value.fdr=pvalues.fdr)
  write.table (data, paste (prefix, '-mrna-cor.tsv', sep=''), row.names=FALSE, sep='\t')
  
  
  ## overall summary and histogram
  s <- summary (correlations)
  s.pos <- sum (correlations > 0, na.rm=TRUE) / length (correlations)
  print (s, quote=FALSE)
  print (table (correlations > 0), quote=FALSE)
  print (paste ('Positive correlations:', round (s.pos*100, digits=2), '%'), quote=F)
  
  pdf (paste (prefix, '-mrna-cor.pdf', sep=''), width=8, height=5)
  q <- qplot (correlations, geom='histogram', binwidth=0.01, xlab=paste (type, '- mrna correlation'))
  print (q)
  dev.off ()
  
  print (table (pvalues.fdr < FDR), quote=FALSE)
  print (paste ('Significant correlations ( FDR <', FDR, '):',
                round (sum (pvalues.fdr < FDR, na.rm=T) / length(pvalues.fdr) * 100, digits=2), '%'), quote=F)
  correlations.sig <- correlations [ pvalues.fdr < FDR ]
  s.sig <- summary (correlations.sig)
  s.sig.pos <- sum (correlations.sig > 0, na.rm=TRUE) / length (correlations.sig)
  print (s.sig, quote=FALSE)
  print (table (correlations.sig > 0), quote=FALSE)
  print (paste ('Positive significant correlations:', round (s.sig.pos*100, digits=2), '%'), quote=F)
  
  pdf (paste (prefix, '-mrna-cor-sig.pdf', sep=''), width=8, height=5)
  q <- qplot (correlations.sig, geom='histogram', binwidth=0.01,
              xlab=paste (type, '- mrna correlation\nFDR pvalue <', FDR))
  print (q)
  dev.off ()
  
  pdf (paste (prefix, '-mrna-cor-combined.pdf', sep=''), width=8, height=5)
  hist (correlations, col='grey', breaks=200, xlab=paste (type, "- mrna correlation"), ylab='count')
  hist (correlations.sig, col='red', breaks=200, add=TRUE)
  legend (-0.5, 80, c ('all', 'significant'), fill=c('grey','red'), bty='n')
  dev.off ()
  
  
  ##
  ## subsets
  ##
  
  
  # profile plots for significant subset
  plot.profiles <- function (data.subset, pdf.file, size) {
    cols.pome <- grep ('\\.pome', colnames (data.subset))
    cols.mrna <- grep ('\\.mrna', colnames (data.subset))
    d.pome <- data.subset [, cols.pome]
    d.mrna <- data.subset [, cols.mrna]
    rownames (d.pome) <- rownames (d.mrna) <- make.names (paste ('Gene:', data.subset[,1], '  Protein:', data.subset[,2]), unique=TRUE)
    
    d.pome <- melt (t (d.pome))
    d.mrna <- melt (t (d.mrna))
    d.plot <- rbind (d.pome, d.mrna)
    colnames (d.plot) <- c ('sample.type', 'gene.gi', 'intensity')
    
    get.sample.name <- function (x) {
      y <- strsplit (toString(x), split='\\.')[[1]]
      paste (y[1:(length(y)-1)], collapse='.')  # to handle replicate with .1 in name
    }
    
    sample <- unlist (lapply (d.plot[,'sample.type'], get.sample.name))
    type <- unlist (lapply (d.plot[,1],
                            function (x) {
                              y <- strsplit (toString(x), split='\\.')[[1]]
                              y [length(y)]
                            }))
    d.plot <- cbind (d.plot, sample, type)
    d.plot[,'sample'] <- factor (d.plot[,'sample'],          # to retain the ordering of samples in plots
                                 levels= unlist (lapply (colnames (data.subset)[cols.pome], get.sample.name)))
    
    q <- ggplot (aes (x=sample, y=intensity, group=type, colour=type), data=d.plot) +
          geom_line() + geom_point() + theme(axis.text.x=element_text (angle=90, vjust=0.5, size=9)) +
          facet_wrap ( ~ gene.gi) 
    
    ggsave (pdf.file, width=size, height=size)
  }
  
  
  # best pairs (high correlation + significant p-value)
  # relax p-value threshold for neg correlation to include some
  select1 <- which ( abs (correlations) > 0.7 & pvalues.fdr < FDR )
  data.subset1 <- data [select1, ]
  write.table (data.subset1, paste (prefix, '-mrna-cor-best.tsv', sep=''), row.names=FALSE, sep='\t')
  plot.profiles (data.subset1, paste (prefix, '-mrna-cor-best.pdf', sep=''), 40)
  
  
  ##
  ## sample-level correlation
  ##
  
  sample.correlations <- sapply (1:length (pome.cols), 
                                 function (i) {
                                   tryCatch (cor.test (data[,pome.cols[i]], data[,mrna.cols[i]], method=correlation.type)$estimate,
                                             error = function (e) NA)
                                 })
  pdf (paste (prefix, '-mrna-sample-cor.pdf', sep=''), width=8, height=5)
  hist (sample.correlations, col='grey', xlab="sample-level protein-mrna correlation", ylab='count')
  dev.off ()
}



## call explore.correlations on the pome with and without bimodal samples
pome.gct.file=file.path (norm.dir, paste (type, '-ratio-norm-nosdfilter-NAmax', na.max, '.gct', sep='')) 
mrna.gct.file='rna-seq.gct' 
pam50.cls=file.path (norm.dir, paste (type, '-pam50.cls', sep=''))
bimodal.cls <- read.cls (file.path (data.dir, bimodal.cls))

explore.correlations (type, pome.gct.file, mrna.gct.file, pam50.cls)

if (sample.source == 'tcga') {
  # whim samples are all unimodal, and unimodal is identical to above
  explore.correlations (paste (type, '-unimodal', sep=''), pome.gct.file, mrna.gct.file, pam50.cls,
                        exclude=(bimodal.cls=='bimodal'))
  ## bimodal samples only
  explore.correlations (paste (type, '-bimodal', sep=''), pome.gct.file, mrna.gct.file, pam50.cls,
                        exclude=(bimodal.cls=='unimodal'))
}



