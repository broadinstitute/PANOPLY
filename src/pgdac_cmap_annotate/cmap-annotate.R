#!/usr/bin/env Rscript

## install pacman and cmapR packages
#  cmapR also requires the prada package which is not automatically installed
if (!require(pacman)) {install.packages("pacman"); library(pacman)}
if (!suppressPackageStartupMessages (require (prada))) install.packages ("prada")
#  make sure 'rhdf5' is loaded BEFORE 'cmapR'
if (!suppressPackageStartupMessages(require(rhdf5))){
  source("https://bioconductor.org/biocLite.R")
  biocLite("rhdf5")
}
if (!suppressPackageStartupMessages(require(cmapR)))
  devtools::install_github("cmap/cmapR")

library (cmapR)
p_load (ggplot2)
p_load (reshape)


## Command line argument processing
## Usage: Rscript cmap-annotate.R tar_file cmap_data_file group data_type enrichment_groups_file output_file config_file
##  tar_file          - has tarball with cmap output results
##  cmap_data_file    - path to annotated_GSE92742_Broad_LINCS_Level5_COMPZ_geneKDsubset_n36720x12328.gctx
##  group             - group for which cmap analysis was run (usually 'all')
##  data_type         - pome, mrna or other data type
##  enrichment_groups_file  - file (similar to experiment design) containing groups for enrichment analysis
##  output_file       - tar file to which output tarball is written
##  config_file       -- (optional) if provided, load parameters from file to over-ride defaults

##  parameter used from config_file
##  alpha             - p-value threshold for cmap profile zscores and enrichments
##  cna_threshold     - copy number up/down if abs (log2(copy number) - 1) is > 0.3
##  log_transfrom     - TRUE/FALSE; should cna data be log transformed?


options( warn = -1 )
args <- commandArgs(trailingOnly=T)

alpha <- 0.05
cna.threshold <- 0.3
log.transform <- FALSE

tar.file <- args[1]
cmap.data.file <- args[2]
group <- args[3]
dtype <- args[4]
groups.file <- args[5]
out.file <- args[6]
if (!is.na (args[7])) source (toString (args[7]))


# extract tarball and determine analysis directory
tar.list <- untar (tar.file, list=TRUE)
analysis.dir <- basename (tar.list[1])
untar (tar.file)

# check to make sure cna and cmap sub-directories exist
if ( !dir.exists (sprintf ("%s/cna", analysis.dir)) || 
     !dir.exists (sprintf ("%s/cmap", analysis.dir)))
  stop ( sprintf ("CNA or CMAP directories not found in tar file %s", tar.file) )


## 
## define functions for annotation and enrichment analysis
##

get.trans.genelist <- function (analysis.dir, cmap.data, grp, typ, alpha) {
  ## extracts trans-genes for enriched CMAP genes and identifies those that are 
  ## extreme in the CMAP profiles and overlap with the trans-genes
  
  cna <- sprintf ("%s/cna", analysis.dir)
  cmap <- sprintf ("%s/cmap", analysis.dir)
  
  ## trans genes for significant genes
  sig.table <- read.delim (sprintf ("%s/%s-cmap-%s-sig-genes-unidirectional.txt", cmap, grp, typ), sep=' ')
  sig.genes <- apply (sig.table, 1, function (x) sprintf ("%s_CNA%s", x[2], tolower (x[1])))
  genesets <- parse.gmt (sprintf ("%s/%s-cmap-%s-updn-genes.gmt", cmap, grp, typ))
  sig.trans <- lapply (sig.genes,
                       function (x) {
                         entry <- as.vector (sapply (genesets[[x]]$entry, function (y) strsplit (y, split=';')[[1]][1]))
                         list (head=x, desc=x, entry=entry)
                       })
  write.gmt (sig.trans, sprintf ("%s/%s-cmap-%s-sig-genes-translist.gmt", cmap, grp, typ))
  
  ## CMAP subsets
  cmap.annot <- read.gctx.meta (cmap.data, dimension='col')
  cmap.subset <- cmap.annot [ cmap.annot[,'pert_iname'] %in% sig.table [,'gene'], ]
  subset.data <- parse.gctx (cmap.data, cid=cmap.subset[,'id'])
  write.gct (subset.data, sprintf ("%s/%s-cmap-%s-sig-genes-profiles.gct", cmap, grp, typ),
             precision=4, appenddim=F, ver=3)
  
  ## CMAP "overlap"
  cmap.nes <- parse.gctx(sprintf ("%s/%s-cmap-%s-nes.gct", cmap, grp, typ))
  outliers <- read.csv (sprintf ("%s/%s-cmap-%s-outliers.csv", cmap, grp, typ))
  overlap <- lapply (1:nrow(sig.table), function (i) {
    print (sig.genes[i])
    o <- outliers [,sig.genes[i]]
    # use make.names below -- there are some ids that have : converted to .
    o.cis <- make.names (cmap.nes@cdesc [cmap.nes@cdesc[,'pert_iname']==sig.table[i,'gene'] & o != 0, 'id'])
    o.data <- data.frame (subset.data@mat [, make.names (subset.data@cid) %in% o.cis])
    o.pval <- data.frame (apply (o.data, 2, function (x) {
      z <- scale (x)
      p <- pnorm (z)
      p < alpha/2 | p > (1-alpha/2)
    }))
    rownames (o.pval) <- subset.data@rid
    o.overlap <- apply (data.frame (o.pval [sig.trans [[i]]$entry, ]), 1, any, na.rm=TRUE)
    list (head=sig.trans[[i]]$head, 
          desc=sprintf ("%d / %d", sum(o.overlap, na.rm=T), length(sig.trans[[i]]$entry)), 
          entry=sig.trans[[i]]$entry[o.overlap])
  })
  write.gmt (overlap, sprintf ("%s/%s-cmap-%s-sig-genes-overlap.gmt", cmap, grp, typ))
  
  ## plot overlap results
  #  trans-genes vs cmap (significant regulated genes) overlap
  d <- overlap
  trans <- as.vector (as.numeric (sapply (d, function (x) strsplit (x$desc, split=" / ")[[1]][2])))
  overlap <- as.vector (as.numeric (sapply (d, function (x) strsplit (x$desc, split=" / ")[[1]][1])))
  trans <- trans - overlap  # since plots are stacked bars
  
  ids <- as.vector (sapply (d, function (x) x$head))
  data <- data.frame (id=factor (ids, levels=ids), trans.gene=trans, cmap.overlap=overlap)
  plot.data <- melt (data)
  colnames (plot.data) <- c ('gene', 'type', 'count')
  p1 <- ggplot (plot.data, aes(y=count, x=gene, fill=type)) + geom_bar(position='stack', stat='identity') + 
    theme(axis.text.x = element_text(angle = 90)) 
  ggsave (sprintf ("%s/%s-cmap-%s-sig-genes-overlap.pdf", cmap, grp, typ))
}



cmap.annot.enrich <- function (analysis.dir,
                               groups.file,
                               group,                # subgroup name on which CNA analysis was run
                               dtype,                # dtype is either pome or mrna
                               cna.threshold,        # copy number up/down if abs (log2(copy number) - 1) is > 0.3
                               log.transform,
                               fdr)
{
  ## identify CNA up/down samples, and calculate enrichment 
  ## for various sample annotations in the groups.file
  
  ## read CNA data, subset to list of CMAP sig genes, and log transform
  cna.dir <- file.path (analysis.dir, 'cna')
  cmap.dir <- file.path (analysis.dir, 'cmap')
  cna <- read.csv (file.path (cna.dir, sprintf ('%s-cna-matrix.csv', group)))
  rownames (cna) <- cna[,1]
  sig.gene.table <- read.delim (sprintf ('%s/%s-cmap-%s-sig-genes-unidirectional.txt', 
                                         cmap.dir, group, dtype), sep=' ')
  sig.genes <- as.character (sig.gene.table [,'gene'])
  cna.data <- cna [sig.genes,2:ncol(cna)]
  if (log.transform) cna.data <- log2 (cna.data) - 1
  # convert to up (+1) / down (-1) calls
  cna.data <- sign (cna.data) * (abs (cna.data) > cna.threshold)
  
  # assemble CNA up/dn sample list
  samples <- colnames (cna.data)
  up <- apply (cna.data [sig.genes [sig.gene.table[,'direction']=='AMP'],], 1, function (x) x==1)
  colnames (up) <- paste (colnames (up), '.AMP', sep='')
  dn <- apply (cna.data[sig.genes [sig.gene.table[,'direction']=='DEL'],], 1, function (x) x==-1)
  colnames (dn) <- paste (colnames (dn), '.DEL', sep='')
  updn <- data.frame (Sample.ID=samples, up, dn)
  write.csv (updn, sprintf ('%s/%s-cmap-%s-sig-genes-extreme-samples.csv', 
                            cmap.dir, group, dtype), row.names=FALSE)
  
  ## plot subset size (ie number of extreme samples)
  s <- updn
  n.extreme <- apply (s[,-1], 2, function (x) sum (x))
  data.ext <- data.frame (gene=factor (names (n.extreme), levels=names (n.extreme)), subset.size=n.extreme)
  p2 <- ggplot (data.ext, aes(x=gene, y=subset.size)) + geom_bar (stat='identity') + 
    theme(axis.text.x = element_text(angle = 90)) 
  ggsave (sprintf ('%s/%s-cmap-%s-sig-genes-extreme-samples.pdf', 
                   cmap.dir, group, dtype))
  
  
  # groups file format similar to expt-design-file with Sample.ID and additional columns
  # enrichement test will be run for each additional column
  subgroup.table <- read.csv (groups.file)
  rownames (subgroup.table) <- subgroup.table [, 'Sample.ID']
  subgroup.table <- subgroup.table [samples,]  # reorder samples (if needed)
  cls.list <- setdiff (colnames (subgroup.table), 'Sample.ID')
  
  if (length (cls.list) > 0) {
    # perform enrichment on sample annotation groups
    enrich <- NULL
    for (cmap.g in colnames (updn)[-1]) {
      cmapg.data <- factor (updn[,cmap.g], c ('TRUE', 'FALSE'))  # TRUE is the class of interest and must be first
      for (g in cls.list) {
        gl <- factor (subgroup.table [,g])
        for (gx in levels (gl)) {
          gx.cls <- ifelse (gl == gx, gx, "other")
          if (any (gx.cls == "", na.rm=T)) next    # skip blank classes
          gx.cls <- factor (gx.cls, levels=c(gx, 'other'))   # class of interest must be first level
          pval <- fisher.test (cmapg.data, gx.cls, alternative='greater')$p.value
          enrich <- rbind (enrich, c (cmap.g, g, gx, pval))
        }
      }
    }
    colnames (enrich) <- c ('cna.gene', 'group', 'subgroup', 'fisher.test.pvalue')
    enrich <- transform (enrich, fisher.test.pvalue=as.numeric (as.character (fisher.test.pvalue)))
    enrich <- cbind (enrich, adj.pvalue=p.adjust (enrich[,'fisher.test.pvalue'], method='BH'))
  }
  write.csv (enrich, sprintf ('%s/%s-cmap-%s-sig-genes-enrichment.csv', 
                              cmap.dir, group, dtype), row.names=FALSE)
  
  # extract statistically significant enrichments
  e <- enrich
  e.sig <- e [ e[,'fisher.test.pvalue'] < fdr, -5]
  write.csv (e.sig, sprintf ('%s/%s-cmap-%s-sig-genes-enrichment-pval%.2f.csv', 
                             cmap.dir, group, dtype, fdr))
}    



create.ssgsea.dataset <- function (analysis.dir, grp, typ) {
  ## create gct dataset for running ssGSEA on trans-genes for CMAP enriched genes
  
  cna <- sprintf ("%s/cna", analysis.dir)
  cmap <- sprintf ("%s/cmap", analysis.dir)
  
  r <- read.csv (sprintf ('%s/%s-%s-vs-cna-corr.csv', cna, grp, typ), row.names=1)

  ## trans p-values for significant genes
  sig.table <- read.delim (sprintf ("%s/%s-cmap-%s-sig-genes-unidirectional.txt", cmap, grp, typ), 
                           sep=' ', stringsAsFactors=FALSE)
  sig.subset <- r [, unique (sig.table[,2])]
  
  # set up GCT v1.3 object
  all.genes <- rownames (sig.subset)
  sig.genes <- colnames (sig.subset)
  file.prefix <- sprintf ('%s/%s-cmap-%s-gsea', cmap, grp, typ)
  gct <- new ('GCT')
  gct@mat <- as.matrix (sig.subset)
  gct@rid <- all.genes
  gct@cid <- sig.genes
  gct@rdesc <- data.frame (GeneSymbol=all.genes, stringsAsFactors=FALSE)
  gct@cdesc <- data.frame (Sample=sig.genes, stringsAsFactors=FALSE)
  gct@version <- "#1.3"
  gct@src <- sprintf ('%s-input.gct', file.prefix)
  # write output file
  write.gct (gct, gct@src, ver=3, precision=4, appenddim=FALSE)
}

# 
#   ## run ssGSEA
#   out.prefix <- sprintf ('%s-output', file.prefix)
#   ssGSEA2 (gct@src, out.prefix, db,
#            sample.norm.type='rank',
#            weight=0,
#            statistic='area.under.RES',
#            output.score.type='NES',
#            nperm=1000,
#            combine.mode='combine.off',
#            min.overlap=10,
#            correl.type='rank',
#            global.fdr=FALSE,
#            extended.output=FALSE,
#            export.signat.gct=FALSE)
#   
#   
#   ## print enrichments
#   p <- parse.gctx ( sprintf ("%s-fdr-pvalues.gct", out.prefix) )
#   e <- apply (p@mat, 2, function (x) p@rid [which (x < FDR)])
#   sink (sprintf ('%s-enrichments.txt', out.prefix))
#   print (e)
#   sink ()
# 


##
## run annotation and enrichment functions
##
get.trans.genelist (analysis.dir, cmap.data.file, group, dtype, alpha)
cmap.annot.enrich (analysis.dir, groups.file, group, dtype, cna.threshold, log.transform, alpha)
create.ssgsea.dataset (analysis.dir, group, dtype)

## create output tar file
tar (out.file, files=analysis.dir)

