if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( cmapR )
p_load( glue )
if ( !requireNamespace("BiocManager", quietly = TRUE ) ) 
  install.packages( "BiocManager" )
BiocManager::install( "org.Hs.eg.db" )
library( org.Hs.eg.db )

main <- function()
{
  csvs <- list.files( path = '.', pattern = "*-analysis-markers-all.csv" )
  
  dir.name <- glue( "{getwd()}/contrasts" )
  if ( !dir.exists( dir.name ) )
    dir.create( dir.name )
  for ( file.name in csvs ) 
  {
    ## Extract contrasts / Fold Change columns as numeric
    csv <- read.csv2( file.name, header = T, stringsAsFactors = F, sep = ',' )
    col <- grep( "contrast*", colnames( csv ), value = T )
    if ( length( col ) == 0 && 'Fold.Change' %in% colnames( csv ) )
      col <- 'Fold.Change'
    mat <- as.data.frame( csv[, col] )
    mat <- as.matrix( apply( mat, c(1, 2), as.numeric ) )
    if ( nrow( mat ) == 0 ) next
    
    ## Keep remainder columns as rdesc
    keep_cols <- setdiff( colnames( csv ), col )
    
    ## Map Protein IDs to gene names
    prot_ids <- csv$Gene.ID
    prot_ids_base <- unlist( lapply( prot_ids, function( x ) 
      unlist( strsplit( x, split = '\\.' ) )[1] ) )
    map.ids <- select( org.Hs.eg.db, keys = prot_ids_base, 
                       columns = c( 'SYMBOL' ), keytype = 'REFSEQ' )
    gene_ids <- unlist( 
      lapply( 1:length( prot_ids ),
              function( x )
              { 
                prot <- prot_ids[x]
                gene <- map.ids$SYMBOL[which( unlist( 
                  strsplit( prot, split = '\\.' ) )[1] == map.ids$REFSEQ )[1]]
                return( gene ) 
              } ) )
    
    ## Format rdesc columns into proper GCT format
    rdesc <- cbind( csv[, keep_cols], geneSymbol = gene_ids )
    rdesc[] <- lapply( rdesc, as.character )
    rownames( mat ) <- csv[, 'Gene.ID']
    colnames( mat ) <- col
    
    ## Annotate column names with more information
    group.name <- unlist( strsplit( file.name, split = '-' ) )[2]
    data.type <- unlist( strsplit( file.name, split = '-' ) )[1]
    cls <- readLines( glue( "{data.type}-{group.name}.cls" ) )[2]
    cls <- unlist( strsplit( cls, split = ' ' ) )[-1]
    col.ann <- c()
    if ( length( col ) > 0 && col[1] != 'Fold.Change' )
    {
      for ( contrast in col )
      {
        cName <- unlist( 
          strsplit( contrast, split = 'contrast.' ) )[2]
        col.ann <- c( col.ann, glue( "{cName}-vs-rest" ) )
      }
    } else if ( col[1] == "Fold.Change" )
      col.ann <- glue( "{cls[2]}-over-{cls[1]}" )
    colnames( mat ) <- col.ann
    print( colnames( mat ) )
    
    target.file.name <- glue( "contrasts/{unlist( strsplit( file.name, split = '-analysis-markers-all.csv' ) )[1]}-contrast.gct" )
    gct_s4object <- new( "GCT",
                         mat = mat,
                         cdesc = data.frame(),
                         rdesc = rdesc,
                         rid = rownames( mat ),
                         cid = colnames( mat ),
                         src = target.file.name )
    write.gct( gct_s4object, target.file.name, ver = 3, appenddim = FALSE )
  }
}

if ( !interactive() ) {
  main()
}

