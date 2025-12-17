preprocessGCT <- function(
    gct.str='',                                           ## path to GCT file 
    level=c('ssc', 'gc', 'gcr'),                          ## single-site-centric, gene-centric or gene-centric-redundant reports
    acc.type=c('uniprot', 'refseq', 'symbol'),            ## accession type of 'rid' object in GCT file
    
    mode=c('mean', 'median', 'sd', 'SGT', 'abs.max'),     ## how should multiple sites per gene be combine; 
    SGT.col='subgroupNum',                                ## column used to collpase to SGT        
    ## sd - most variable (standard deviation) across sample columns
    ## SGT - subgroup top: first subgroup in protein group (spectrum mill)
    ## abs.max - most extreme site; for log-transformed, signed p-values
    
    ## gene centric parameters
    gene.col='geneSymbol',                                ## name of column listing gene names; used for gene centric reports and to fix gene names 
    fix.gene.names=F,                                     ## if TRUE, try to fix the gene2date conversion problem; works for site-centric reports as well
    fix.gene.names.column='id.description',               ## alternative column in rdesc containing gene names that
    ## can be used to fix gene names, e.g. description: 'septin-10 isoform 4 GN=SEPT10'
    fix.gene.names.regexpr='.* GN\\=(.*)$',               ## regular expression to extract gene names from 'fix.gene.names.column'
    humanize.gene.names=FALSE,                            ## if TRUE, gene symbols will be capitalized (for e.g. mouse or rat)
    
    ## ssc specific parameters
    id.type=c('sm', 'wg', 'ph'),                          ## notation of site-ids: sm-Spectrum Mill; wg-Web Gestalt; ph-Philosopher
    id.type.out=c('uniprot', 'refseq', 'seqwin', 'psp'),  ## type of site id for output; psp not implemented yet
    loc=T,                                                ## if TRUE only fully localized sites will be considered
    seqwin.col='VMsiteFlanks',                            ## column containing flanking sequences, separated by '|'. Only relevant
    mod.res=c('S|T|Y', 'K'),                              ## modified residue(s) 
    mod.type=c('p', 'ac', 'ub'),                          ## modification type
    
    ## functional params
    appenddim=T,                                          ## see cmapR::write.gct()
    preprocess.gct=T                                      ## flag, if FALSE nothing will be done; probably needed for to make this step optional in a FireCLoud WDL.
    
) {
  
  ## immediately return
  if(!preprocess.gct){
    return(gct.str)
  }
  
  require(cmapR)
  require(magrittr)
  require(glue)
  require(dplyr)
  
  ## Process Parameters ---------------
  
  level <- match.arg(level)
  mode <- match.arg(mode)
  mod.res <- match.arg(mod.res)
  mod.type <- match.arg( mod.type )
  acc.type <- match.arg( acc.type )
  #org <- match.arg(org)
  id.type <- match.arg(id.type)
  id.type.out<- match.arg(id.type.out)
  
  ## Import GCT ---------------
  
  if(file.exists(gct.str)){
    cat('importing gct file: ', gct.str, ' ...\n')
    gct <- try(parse.gctx(gct.str))
  } else {
    stop(glue("File '{gct.str}' not found!\n"))
  }
  if(class(gct) == 'try-error'){
    
    ## - cmapR functions stop if ids are not unique
    ## - import gct using readLines and make ids unique
    if(length(grep('rid must be unique', gct) ) > 0) {
      gct.tmp <- readLines(gct.str)
      #first column
      rid <- gct.tmp %>% sub('\t.*','', .)
      #data and meta data columns
      meta <- strsplit(gct.tmp[2], '\t') %>% unlist() %>% as.numeric()
      rid.idx <- (meta[4]+3) : length(rid)
      #check whether ids are unique
      if(length(rid[rid.idx]) > length(unique(rid[rid.idx]))){
        warning('rids not unique! Making ids unique and exporting new GCT file...\n\n')
        #make unique
        rid[rid.idx] <- make.unique(rid[rid.idx])
        #other columns
        rest <- gct.tmp %>% sub('.*?\t','', .)
        rest[1] <- ''
        gct.tmp2 <- paste(rid, rest, sep='\t') 
        gct.tmp2[1] <-  sub('\t.*','',gct.tmp2[1])
        #export
        gct.unique <- sub('\\.gct', '_unique.gct', gct.str)
        writeLines(gct.tmp2, con=gct.unique)
        
        gct <- parse.gctx(fname = gct.unique)
      }
    } #end if 'rid not unique'
  }
  
  ## Preprocess IDs ---------------
  if (level %in% c('gc', 'gcr')) {
    ### Gene-Centric Preprocessing ---------------
    #### fix gene names affected by gene2date conversion (optional) ----
    if(fix.gene.names){
      
      if(!gene.col %in% colnames(gct@rdesc))
        stop(glue("Column {gene.col} not found!"))
      
      if(!fix.gene.names.column %in% colnames(gct@rdesc))
        stop(glue("Column {fix.gene.names.column} not found!"))
      
      genes <- sub(fix.gene.names.regexpr, '\\1', gct@rdesc[, fix.gene.names.column])
      gct@rdesc[, gene.col] <- genes  
    }
    
    #### humanize gene symbols (optional) ----
    if(humanize.gene.names){
      genes <- toupper(genes)
    }
    
    #### remove empty gene symbols ----
    keep.idx <- !(is.na(gct@rdesc[, gene.col]) | nchar(gct@rdesc[, gene.col]) == 0) # logical vector of rows without NA or empty genesymbols
    gct.filt = subset_gct(gct, rid = keep.idx)
    if (sum(keep.idx) != length(keep.idx))
      warning(glue("Removed {sum(!keep.idx)} rows due to missing gene symbol in column '{gene.col}'\n\n"))
    
    #### GCR / GC specific wrangling ----
    if (level == "gcr") {
      ##### GCR Wrangling ----
      gct.fin = gct.filt
      gct.fin@rdesc$id.original = gct.fin@rid # retain original ID
      gct.fin@rid = make.unique(gct.fin@rdesc[, gene.col], sep = '_') # overwrite rid with (unique) gene symbols
    } else {
      ##### GC Collapsing ----
      
      ## aggregate data
      if(mode == 'mean') {
        ###### mean -----
        mat.gc <- aggregate(gct.filt@mat, FUN=function(x) mean(x, na.rm=T), by=list(gct.filt@rdesc[, gene.col]))  
        rdesc.tmp = NULL # aggregate rdesc using default method
      } else if(mode == 'median'){
        ###### median -----
        mat.gc <- aggregate(gct.filt@mat, FUN=function(x) median(x, na.rm=T), by=list(gct.filt@rdesc[, gene.col]))  
        rdesc.tmp = NULL # aggregate rdesc using default method
      } else if(mode == 'SGT'){
        ###### SGT -----
        
        # 1. extract top subgroup
        # 2. collapse to genes by taking median expression
        if(!SGT.col %in% colnames(gct.filt@rdesc)) stop(glue("Column {SGT.col} not found!"))
        
        ## keep lowest subgroup per gene
        sgt.idx <- tapply(gct.filt@rdesc[, SGT.col], # SGT column
                          gct.filt@rdesc[, gene.col], # geneSymbol column
                          function(x) x[ which.min( as.numeric(sub('.*\\.', '', x))) ] )
        
        ## remove duplicated subgroup numbers introduced when there are more than 9 subgroups (happens in parsing module of panoply)
        ## 234.10 -> 234.1
        sgt.idx <- unique(sgt.idx)
        
        keep.idx <- match(sgt.idx, gct.filt@rdesc[, SGT.col])
        
        mat.gc <- data.frame(id=gct.filt@rdesc[, gene.col][keep.idx], # first column with gene symbols
                             gct.filt@mat[keep.idx,,drop=F]) # data to keep
        rdesc.tmp = gct.filt@rdesc[keep.idx,,drop=F] # keep relevant rdesc values
      } else  if(mode == 'abs.max'){
        ###### abs.max -----
        mat.gc <- aggregate(gct.filt@mat, FUN=function(x)x[ which.max(abs(x)) ], by=list(gct.filt@rdesc[, gene.col]))   
        rdesc.tmp = NULL # aggregate rdesc using default method
      }
      
      ##### aggregate GCT object -----
      rid <- mat.gc[, 1] %>% as.character
      mat.gc <- mat.gc[, -c(1), drop=F]
      ## collect rdesc
      if (!is.null(rdesc.tmp)) { # custom rdesc aggregation
        rdesc.gc = rdesc.tmp 
      } else { # standard rdesc aggregation
        rdesc.gc <- aggregate(gct.filt@rdesc, FUN=function(x) paste(unique(x), collapse='|'), by=list(gct.filt@rdesc[, gene.col]))
        rdesc.gc <- rdesc.gc[, -c(1), drop=F]
      }
      ## apply rid as rownames
      rownames(mat.gc) <- rownames(rdesc.gc) <- rid
      ## create GCT object
      gct.fin = GCT(mat = data.matrix(mat.gc),
                    rid = rid,
                    cid = colnames(mat.gc),
                    rdesc = rdesc.gc,
                    cdesc = gct.filt@cdesc)
      
      # sanity check finalized GCT file
      if ( any(gct.fin@rid!= rownames(gct.fin@mat)) ) stop('Mismatch found between rid and matrix row-names in final GCT')
      if ( any(gct.fin@rid!= rownames(gct.fin@rdesc)) ) stop('Mismatch found between rid and rdesc rownames in final GCT')
      if ( any(gct.fin@cid!= colnames(gct.fin@mat)) ) stop('Mismatch found between cid and matrix column-names in final GCT')
      if ( any(gct.fin@cid!= rownames(gct.fin@cdesc)) ) stop('Mismatch found between cid and cdesc rownames in final GCT')
      
    }
    
  } else if (level == 'ssc') {
    ### Site-Centric Preprocessing ---------------
    
    
    #### Parse Site Information ---------------
    
    
    #### Localization ---------------
    
    
    #### Multi-Site Processing ---------------
    
    
    
    #### Aggregate Duplicate Sites ---------------
    
  } else {
    stop(glue("Provided level '{level}' is not a valid option."))
  }
  
  ## Create Output File ---------------
  
}
    