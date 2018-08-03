
library (cmapR)

## create subset of CMAP dataset with only gene knockdown profiles
## split the dataset with gene knockdown profiles into subsets with about 400-500 profiles each

input.cmap.ds <- 'annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx'   # CMAP Touchstone dataset
output.subset.ds <- 'annotated_GSE92742_Broad_LINCS_Level5_COMPZ_geneKDsubset'
subsets.dir <- 'cmap-data-subsets'
N <- 100    # number of subsets


# read column annotations in the Touchstone dataset, and extract gene knockdown profiles
col.annot <- read.gctx.meta (input.cmap.ds, dimension = 'col')
gene.kds <- col.annot [ col.annot[,'pert_type']=="trt_sh.cgs", ]

# extract and save gene knockdown profiles
# ensure that row ids are gene symbols (required for ssGSEA)
kd.subset <- gene.kds [,'id']
kd.data <- parse.gctx (input.cmap.ds, cid = kd.subset)
kd.data@rid <- kd.data@rdesc[,'pr_gene_symbol']
kd.data@rdesc[,'id'] <- kd.data@rid

write.gctx (kd.data, output.subset.ds)


# split into smaller subsets and save into subsets.dir
# create an index file to track subsets (for WDL scatter)
if (! dir.exists (subsets.dir)) dir.create (subsets.dir)
subset.size <- ceiling (ncol (kd.data@mat) / N)
index <- NULL
for (i in 1:N) {
  index.start <- (i-1) * subset.size + 1
  index.end <- min (i * subset.size, ncol(kd.data@mat))

  subset.i <- subset.gct (kd.data, cid=index.start:index.end)
  subset.file.name <-  sprintf ("%s-%03d.gctx", output.subset.ds, i)
  write.gctx (subset.i, file.path (subsets.dir, subset.file.name), appenddim=FALSE)
  
  if (index.end >= ncol(kd.data@mat)) break
  else index <- c (index, subset.file.name)
}
write.table (index, 'cmap-data-subsets-index.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)


## also write out list of genes knocked down in the CMAP data
#  this is used in selecting query genes when creating CMAP input
cmap.kd.genes <- unique (gene.kds [,'pert_iname'])
write.table (cmap.kd.genes, 'cmap-knockdown-genes-list.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)

