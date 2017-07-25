
# $Id$
print ("$Id")


source ('preamble.r')
Source ('markersel-classify.r')
min.mutation.N <- 5


write.R.file <- function (gene, gct, cls) {
  # create R script to run marker selection with
  # (i) FDR=0.01, GSEA enabled (defaults)
  # (ii) FRD=0.1, GSEA disabled
  script <- paste ('gene-', gene, '.r', sep='')
  run.gsea <- ifelse (type=='proteome', TRUE, FALSE)   # can't run GSEA for phosphotproteome unless gene mapping is finalized
  cat ("source ('preamble.r')\n",
       "Source ('markersel-classify.r')\n",
       "marker.selection.and.classification ('", gct, "', '", cls, 
                                             "', 'gene-", gene, "/", gene, "', gsea.prefix='", gene,
                                             "', id.col=", id.col, ", desc.col=", desc.col, ", gsea=", run.gsea,
                                             ", duplicate.gene.policy='", duplicate.gene.policy,
                                             "', output.imputed.data=TRUE, output.all.markers=TRUE)\n",
       "marker.selection.and.classification ('", gct, "', '", cls, 
                                             "', 'gene-", gene, "/", gene, "', gsea.prefix='", gene,
                                             "', id.col=", id.col, ", desc.col=", desc.col, 
                                             ", fdr=0.1, gsea=FALSE, output.imputed.data=FALSE, output.all.markers=FALSE)\n",
       file=script, sep='')
}


run.R.file <- function (gene) {
  script <- paste ('gene-', gene, '.r', sep='')  
  cmd <- paste ('bsub -q priority -o lsf.out -R "rusage[mem=8000]" R CMD BATCH --vanilla', script)
  system (cmd)
}


## read and filter table of significant / major mutated genes
d <- read.csv (file.path (data.dir, 'major-mutated-genes.csv'), check.names=FALSE)

mut.table <- d[, grep (sample.source, colnames(d), ignore.case=TRUE)]
rownames (mut.table) <- d [, 'GeneSymbol']
if (sample.source == 'tcga') {
  colnames (mut.table) <- sapply (colnames (mut.table), 
                                  function (x) strsplit (x, split='TCGA-')[[1]][2])
}
    

## read in p-ome data
pome <- read.gct (master.file, check.names=FALSE)
pome.data <- pome [, 3:ncol(pome)]
samples <- sapply (colnames (pome.data), function (x) strsplit (x, split='\\.')[[1]][1])

## eliminate items in p-ome data not present in mut.table
present <- sapply (samples, function (x) any (grepl (x, colnames (mut.table))))
pome.data <- pome.data [, present]
samples <- samples [ present ]
order <- sapply (samples, function (x) grep (x, colnames (mut.table)))
if (sample.source == 'whim') {
  # whim samples need to be matched using the whole name
  order <- sapply (samples, function (x) grep (paste ('^', x, '$', sep=''), colnames (mut.table)))
}


## write out new gct file
new.gct <- file.path (norm.dir, paste (master.prefix, '-brca-samples.gct', sep=''))
write.gct (cbind (pome[,1:2], pome.data), new.gct)


## for each gene in the mutation table, create cls file and run marker selection
for (i in 1:nrow(mut.table)) {
  gene <- rownames (mut.table)[i]
  if ( sum (mut.table[i, order]) >= min.mutation.N) {
    # (if there are enough mutated samples)
    new.cls <- paste ('gene-', gene, '.cls', sep='')
    cls <- sapply (mut.table [i, order], ifelse, "mutated", "unmutated")
    write.cls (cls, new.cls)
    
    gene.dir <- paste ('gene-', gene, sep='')
    if (!file.exists (gene.dir)) dir.create (gene.dir)
    write.R.file (gene, new.gct, new.cls)
    run.R.file (gene)
  } else {
    cat ('Not enough mutated samples: Skipping analysis for gene', gene, '... \n')
  }
}
