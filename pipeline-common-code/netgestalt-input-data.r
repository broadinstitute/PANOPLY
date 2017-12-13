

source ('preamble.r')
Source ('misc.r')
Source ('map-to-genes.r')


# indicate what data to process
# (for selective processing of specific data)
# process.all <- FALSE      # overrides individual settings
process.all <- TRUE      # overrides individual settings
process.rnaseq <- TRUE
process.pome <- TRUE
process.cna <- TRUE
process.association <- TRUE


# gene symbol is in the "Description" column of gct files (RNA, proteome)
# 'Name' column for CNA
rna.gene.id.col <- pome.gene.id.col <- 'Description'
cna.gene.id.col <- 'Name'


##
## Samples
##
# get list of samples in the correct order -- this order will be used for all datasets
# samples <- scan (file.path (norm.dir, paste (master.prefix, 'gct', sep='.')),
#                  what='char', nlines=1, skip=2)[c(-1,-2)]

# # define bimodal / normal samples
# bimodality.coeff.cls <- read.cls (file.path (data.dir, bimodal.cls))
# 
# # determine whether to process all samples or just unimodal samples
# if (data.subset == "unimodal") {
#   samples <- samples [ bimodality.coeff.cls == "unimodal" ]
# }

 
data.matrix <- function (gct.file, gene.id.col) {
  d <- read.gct2 (gct.file)
  s <- colnames (d) [c(-1,-2)]
  m <- d [, c (gene.id.col, s)]
  colnames (m) <- c ('GeneSymbol', s)
  d.matrix <- process.duplicate.genes (m, genesym.col=1, data.cols=2:ncol(m),
                                       official.symbols.file=file.path (data.dir, official.genesyms),
                                       map.genes=TRUE)
  return (list (matrix=d.matrix, samples=s))
}


##
## RNA-seq Data
##
if (process.all || process.rnaseq) {
  rna <- data.matrix (file.path (rnaseq.dir, 'rna-seq.gct'), rna.gene.id.col)
  # m <- read.gct2 ( file.path (rnaseq.dir, 'rna-seq.gct'), check.names=FALSE )
  # rna.samples <- colnames (m) [c(-1,-2)]
  # rnaseq.matrix <- m [, c (gct.gene.id.col, rna.samples)]
  # colnames (rnaseq.matrix) <- c ('GeneSymbol', samples)
  # rnaseq.matrix <- process.duplicate.genes (rnaseq.matrix, genesym.col=1, data.cols=2:ncol(rnaseq.matrix),
  #                                           official.symbols.file=file.path (data.dir, official.genesyms),
  #                                           map.genes=TRUE)
  # write.csv (rnaseq.matrix, 'rnaseq-matrix.csv', row.names=FALSE)
}




##
## Global Proteome or Phosphoproteome Data
##
if (process.all || process.pome) {
  pome.dir <- norm.dir
  pome <- data.matrix (file.path (pome.dir, paste (master.prefix, '.gct', sep='')), pome.gene.id.col)
  # pr.n <- read.gct2 ( file.path (prot.dir, paste (master.prefix, '.gct', sep='')), check.names=FALSE )
  # pr <- cbind (pr.n[,gct.gene.id.col], pr.n[,3:ncol(pr.n)])
  # 
  # proteome.matrix <- pr
  # colnames (proteome.matrix)[1] <- 'GeneSymbol'
  # proteome.matrix <- process.duplicate.genes (proteome.matrix, genesym.col=1, data.cols=2:ncol(proteome.matrix),
  #                                             official.symbols.file=file.path (data.dir, official.genesyms),
  #                                             map.genes=TRUE, policy=duplicate.gene.policy)
  # write.csv (proteome.matrix, paste (type, '-matrix.csv', sep=''), row.names=FALSE)
}

  

##
## CNA Data
##
if (process.all || process.cna) {
  cna <- data.matrix (file.path (cna.dir, cna.data.file), cna.gene.id.col)
  #   # CNA data in a single file
  #   cna.data <- read.gct2 (file.path (cna.dir, cna.data.file))
  #   cna.data <- cbind (cna.data[,cna.gene.id.col], cna.data[,3:ncol(cna.data)])
  #   cna.samples <- sapply (colnames (cna.data)[-1], 
  #                          function (x) paste ('MB',
  #                                              sprintf ("%03d",
  #                                                       as.numeric (strsplit (sub ('BL', '', x), split='MB')[[1]][2])),
  #                                              sep=''))
  #   colnames (cna.data) <- c ('GeneSymbol', cna.samples)
  #   
  #   # reorder samples to be in the order in which they were run
  #   cna.matrix <- cna.data [, c ('GeneSymbol', samples)]
  # 
  # 
  # # CNA data has gene symbols from multiple sources (and too many genes ...)
  # # keep only those genes that are present either in mrna/rna-seq or proteome
  # #  cna.matrix <- map.genesym (cna.matrix)   
  # # this is too inefficient -- keep only official symbols
  # cna.genes <- toupper (as.character (cna.matrix [,'GeneSymbol']))
  # if (process.all) {
  #   keep <- cna.genes %in% mrna.matrix[,'GeneSymbol'] | cna.genes %in% rnaseq.matrix[,'GeneSymbol'] | cna.genes %in% proteome.matrix[,'GeneSymbol']
  #   cna.matrix <- cna.matrix [ keep, ]
  #   cna.matrix <- process.duplicate.genes (cna.matrix, genesym.col=1, data.cols=2:ncol(cna.matrix),
  #                                          map.genes=FALSE)
  # } else {
  #   cna.matrix <- process.duplicate.genes (cna.matrix, genesym.col=1, data.cols=2:ncol(cna.matrix),
  #                                          official.symbols.file=file.path (data.dir, official.genesyms),
  #                                          map.genes=TRUE)
  # }
  # 
  # write.csv (cna.matrix, 'cna-matrix.csv', row.names=FALSE)
}



##
## Mutation Status Data
##  (not available for WHIM samples)
##
# if (sample.source != "whim" && (process.all || process.mutation) ) {
#   mut.matrix <- read.csv ( file.path (mut.dir, 'mutated-genes-binary.csv'), check.names=FALSE )[,-1]
#   colnames (mut.matrix)[1] <- 'GeneSymbol'
#   
#   # map gene symbols to standard gene names
#   mut.matrix <- process.duplicate.genes (mut.matrix, genesym.col=1, data.cols=2:ncol(mut.matrix),
#                                          official.symbols.file=file.path (data.dir, official.genesyms),
#                                          policy="union")
#   write.csv (mut.matrix, 'mutation-matrix.csv', row.names=FALSE, quote=FALSE)
# }


##
## Derive common sample subset and write out files
##
samples <- NULL
if (exists ("rna")) samples <- rna$samples
if (exists ("pome")) {
  if (is.null (samples)) samples <- pome$samples
  else samples <- intersect (samples, pome$samples)
}
if (exists ("cna")) {
  if (is.null (samples)) samples <- cna$samples
  else samples <- intersect (samples, cna$samples)
}


if (!is.null (samples)) {
  
  write.subset <- function (d, f) {
    d.matrix <- d[, c ('GeneSymbol', samples)]
    write.csv (d.matrix, f, row.names=FALSE, quote=FALSE)
  }
  
  if (exists ("rna")) write.subset (rna$matrix, "rnaseq-matrix.csv")
  if (exists ("pome")) write.subset (pome$matrix, paste (type, '-matrix.csv', sep=''))
  if (exists ("cna")) write.subset (cna$matrix, "cna-matrix.csv")
  
  
  ##
  ## Sample Info
  ##
  sinfo <- read.csv ( file.path (data.dir, 'exptdesign.csv') )
  rownames (sinfo) <- sinfo [,'Sample.ID']
  sinfo <- sinfo [ samples, ]

  # write out properly formatted info file
  out.file <- 'sample-info.csv'
  colnames (sinfo) [1] <- 'Barcode'
  col.names <- colnames (sinfo)
  col.types <- c ('data_type', rep ('CAT', ncol(sinfo)-1))
  write.csv (col.names, out.file, row.names=FALSE, col.names=FALSE)
  write.csv (col.types, out.file, append=TRUE, row.names=FALSE, col.names=FALSE)
  write.csv (sinfo, out.file, append=TRUE, row.names=FALSE, col.names=FALSE)
  
  # recreate cls files (since sample subset included in matrix files could be different)
  for (i in 2:ncol(sinfo)) {
    write.cls (sinfo[,i], sprintf ("%s.cls", colnames (sinfo)[i]))
  }
}
