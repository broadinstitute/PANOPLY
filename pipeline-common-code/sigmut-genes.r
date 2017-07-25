
source ('preamble.r')


# get list of significantly mutated genes (from TCGA breast cancer paper)
tcga.variants <- read.csv ( file.path (data.dir, 'tcga-sigmut-genes.csv') )
tcga.variants[,'GeneSymbol'] <- as.character (tcga.variants[,'GeneSymbol'])
genes <- tcga.variants [,'GeneSymbol']


# extract these genes from the proteome table
data <- read.csv ( file.path (netgestalt.dir, paste (type, '-matrix.csv', sep='')), check.names=FALSE)
data[,'GeneSymbol'] <- as.character (data[,'GeneSymbol'])
keep <- sapply (data[,'GeneSymbol'], function (x) x %in% genes)
data.subset <- data [keep, ]
rownames (data.subset) <- data.subset [,'GeneSymbol']

genes.subset <- merge (tcga.variants[,c(1,4)], data.frame (GeneSymbol=data.subset [,'GeneSymbol']))
data.subset <- data.subset [ genes.subset [order (genes.subset[,2], decreasing=TRUE), 'GeneSymbol'],]


# rearrange columns to be in iTRAQ run order
cnames <- colnames (read.gct (master.file, check.names=FALSE))[c(-1,-2)]
samples <- sapply (cnames, function (x) strsplit (x, split='\\.')[[1]][1])


# retain only BC samples
cls <- read.cls (file.path (norm.dir, paste (type, subset.str, '-pam50.cls', sep='')))
keep <- cls == "Basal" | cls == "Her2" | cls == "LumA" | cls == "LumB"
samples <- samples [keep]   


# write out final table
data.subset <- data.subset [, c('GeneSymbol', samples)]
write.csv (data.subset, 'sigmut-genes.csv', row.names=FALSE, quote=FALSE)

