#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

source ('config.r')

if (!require("pacman")) install.packages ("pacman")
pacman::p_load (RColorBrewer, ComplexHeatmap, circlize, cluster, ape, reshape, ggplot2)
if (! require (estimate)) {
  # only available at R forge
  install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)
  library (estimate)
}


rna.dir <- '../rna'
cor.threshold <- 0.4
heatmap.colors <- brewer.pal (11, "RdYlBu")[2:10]
samplecolors <- c ( brewer.pal (8, name='Dark2'), brewer.pal (9, name='Set1')[-6] )


read <- function (f, suffix) {
  d <- read.csv (f, row.names=1)
  colnames (d) <- paste (colnames(d), suffix, sep='.')
  return (d)
}


## read data and filter by protein-RNA correlation
rna <- read ( sprintf ('%s/rna-matrix.csv', harmonize.dir), 'RNA' )
pome <- read ( sprintf ('%s/%s-matrix.csv', harmonize.dir, type), 'PROT')
cna <- read ( sprintf ('%s/cna-matrix.csv', harmonize.dir), 'CNA')

cis.cor.file <- sprintf ('%s/proteome-mrna-cor.tsv', rna.dir)
if (file.exists(cis.cor.file)) {
  cis.cor <- read.delim (cis.cor.file)
  rownames (cis.cor) <- cis.cor [,gene.id.col] #replaced 'geneSymbol' with gene.id.col
  genes <- intersect (intersect (rownames (cis.cor), rownames (pome)), intersect (rownames (cna), rownames(rna)))
  
  cis.cor <- cis.cor [ genes, ]
  keep <- cis.cor [,'correlation'] > cor.threshold & !is.na (cis.cor[,'correlation']) 
  keep.genes <- cis.cor [keep ,gene.id.col] #replaced 'geneSymbol' with gene.id.col
  
  cna.cx <- cna [keep, ]
  rna.cx <- rna [keep, ]
  pome.cx <- pome [keep, ]
} else {
  warning ('No RNA-protein correlation found -- using ALL genes/proteins')
  cna.cx <- cna
  rna.cx <- rna
  pome.cx <- pome
}


pdf ('sample-qc-plots.pdf', width=12, height=10, pointsize=8)

## correlations (heatmap)
draw.heatmap <- function (d, title) {
  heatmap.range <- range (as.vector (unlist (d)), na.rm=TRUE)
  h <- Heatmap (d, col=colorRamp2 (c(heatmap.range[1], 0, heatmap.range[2]), heatmap.colors [c(11, 6, 1)]),
                cluster_columns=FALSE, cluster_rows=FALSE, name=title,
                heatmap_legend_param=list (color_bar='continuous'), 
                row_names_gp = gpar (fontsize=6), column_names_gp = gpar(fontsize=6))
  draw(h)
}

pVr <- cor (pome.cx, rna.cx, use="pairwise.complete.obs", method="spearman")
draw.heatmap (pVr, "Proteome vs RNA Correlation")
pVc <- cor (pome.cx, cna.cx, use="pairwise.complete.obs", method="spearman")
draw.heatmap (pVc, "Proteome vs CNA Correlation")
rVc <- cor (rna.cx, cna.cx, use="pairwise.complete.obs", method="spearman")
draw.heatmap (rVc, "RNA vs CNA Correlation")


## fan plot
#  (for protein and RNA only -- CNA data is too different)
draw.fanplot <- function (d, title, n.types=2) {
  dist.spearman <- 1 - cor (d, method="spearman", use="pairwise.complete")
  joint.cluster <- agnes (dist.spearman, diss=TRUE, method='complete')
  joint.phy <- as.phylo (as.hclust (joint.cluster))
  # graphics parameters
  n <- ncol (d) / 2
  xmax <- max (joint.phy$edge.length, na.rm=TRUE) * 6 / n.types
  cmax <- n / ceiling (n / length (samplecolors))
  # create plot
  plot (joint.phy, type='fan', tip.color = samplecolors[1:cmax], no.margin=FALSE, cex=0.7, main=title, x.lim=c(-xmax, xmax))
  if (n.types != 1)  ring (0.025, joint.phy, offset=0.2,  col=rep (c ('grey80', 'grey40'), each=n))
}

pAr <- data.frame (pome.cx, rna.cx)
draw.fanplot (pAr, "Proteome vs RNA", n.types=2)


## estimate scores
runEstimate <- function (ds, out.file, type) {
  # run estimate and return scores
  tmp.file <- sprintf ("tmp-%s.gct", Sys.getpid())
  write.gct2 (data.frame (Name=rownames (ds), Description=rownames(ds), ds), tmp.file)
  estimateScore (tmp.file, out.file)
  file.remove (tmp.file)
  
  scores.table <- parse.gctx (out.file)
  scores <- t (scores.table@mat)
  #colnames (scores) <- paste (colnames (scores), type, sep='')
  result <- data.frame (sample=sapply (rownames(scores),
                                       function (x) strsplit (x, split=type)[[1]][1]), type=type, scores)
  
  return (result)
}

rna.ES <- runEstimate (rna.cx, 'rna-estimate-scores.gct', '.RNA')
cna.ES <- runEstimate (cna.cx, 'cna-estimate-scores.gct', '.CNA')
pome.ES <- runEstimate (pome.cx, 'pome-estimate-scores.gct', '.PROT')

results.ES <- rbind (rna.ES,cna.ES,pome.ES)
plot.data <- melt (results.ES, id.vars=c('sample','type'))
ggplot (aes (x=type, y=value, group=type, color=type), data=plot.data) +
  geom_boxplot() + facet_wrap (~ variable, scales='free') +
  ggtitle ('ESTIMATE scores for RNA, CNA and Protein')


dev.off()



