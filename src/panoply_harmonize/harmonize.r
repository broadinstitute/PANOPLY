#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#


source ('config.r')
Source ('misc.r')
Source ('gct-io.r')
Source ('map-to-genes.r')
Source ('process_duplicate_genes.r')

##
## harmonized RNA, CNA and PTM/proteome data by selecting common
## rows (genes) and columns (samples)
## NB: All datasets should have the same samples names (order can be different)
## RNA and CNA should have gene symbols as id's (Name column in gct v1.2/3)
## Parsed PTM/proteome (see filter.r) will have gene symbols in Description column (gct v1.2/3)
##
## output tables are txt format and will have only a GeneSymbol column in addition to the data
## these tables can be used as input to NetGestalt of other gene-centric analysis tool
##

data.matrix.harmonize <- function( gct.file, gene.id.col, policy ){
  # read in data and eliminate rows ("genes") with too many missing or 0 values
  d <- parse.gctx( gct.file )
  remove <- apply( d@mat, 1, function( x ){
    n <- length( x )
    zero.frac <-  sum( x == 0, na.rm = TRUE )/n
    na.frac <- sum( is.na( x ) )/n
    return ( zero.frac > na.max || na.frac > na.max ) } )
  d <- row.subset.gct( d, !remove )
  s <- d@cid
  if( policy == "SGT" ){
    m <- data.frame( GeneSymbol=d@rdesc[, gene.id.col],
                     subgroupNum=d@rdesc[, "subgroupNum"],
                     d@mat )
    data.cols.start <- 3 }
  else{
    m <- data.frame( GeneSymbol = d@rdesc[, gene.id.col], d@mat )
    data.cols.start <- 2 }
  d.matrix <- process.duplicate.genes.2(
    m,
    genesym.col = 1,
    data.cols = data.cols.start:ncol( m ),
    official.symbols.file = file.path( data.dir, official.genesyms ),
    map.genes = TRUE,
    policy = policy )
  rownames( d.matrix ) <- d.matrix [,'GeneSymbol']
  return( list( matrix = d.matrix, samples = s, col.descr = d@cdesc ) )
}


## read in data tables
pome.gene.id.col <- 'GeneSymbol'
rna.gene.id.col <- cna.gene.id.col <- 'id'
# RNA-seq Data
rna <- data.matrix.harmonize( rna.data.file, rna.gene.id.col,
                              duplicate.gene.policy )
# PTM/Proteome Data
if( !( exists( "pome.duplicate.gene.policy" ) ) )
  pome.duplicate.gene.policy <- duplicate.gene.policy
pome <- data.matrix.harmonize( file.path( norm.dir, master.file ),
                               pome.gene.id.col,
                               pome.duplicate.gene.policy )
# CNA Data
cna <- data.matrix.harmonize( cna.data.file, cna.gene.id.col,
                              duplicate.gene.policy )



## harmonize
common.samples <- intersect (pome$samples, intersect (rna$samples, cna$samples))
common.genes <- intersect (rownames (pome$matrix), intersect (rownames (rna$matrix), rownames (cna$matrix)))

if (length (common.samples) > 0 && length (common.genes) > 0) {

  write.subset <- function (d, f) {
    d.matrix <- d[common.genes, c ('GeneSymbol', common.samples)]
    write.csv (d.matrix, f, row.names=FALSE, quote=FALSE)
  }

  write.subset (rna$matrix, "rna-matrix.csv")
  write.subset (pome$matrix, paste (type, '-matrix.csv', sep=''))
  write.subset (cna$matrix, "cna-matrix.csv")

  ## Sample Info
  #  obtain from exptdesign file or from column descriptors in gct v1.3 pome
  if (file.exists (expt.design.file) || nrow (pome$col.descr) > 0) {
    if (file.exists (expt.design.file)) {
      sinfo <- read.csv (expt.design.file)
    } else {
      sinfo <- pome$col.descr
    }
    rownames (sinfo) <- sinfo [,'Sample.ID']
    sinfo <- sinfo [ common.samples, ]
    if (any (is.na (sinfo [,'Sample.ID']))) stop ("Sample information missing for some samples -- check experiment design file")

    # write out properly formatted info file
    out.file <- 'sample-info.csv'
    write.csv (sinfo, out.file, row.names=FALSE, quote=FALSE)

    # recreate cls files (since sample subset included in matrix files could be different)
    for (g in setdiff (colnames (sinfo), c('id', 'Sample.ID'))) {
      write.cls (sinfo[,g], sprintf ("%s.cls", g))
    }
  }
} else {
  stop ("Harmonized datasets have no data -- check data formats and sample names")
}

