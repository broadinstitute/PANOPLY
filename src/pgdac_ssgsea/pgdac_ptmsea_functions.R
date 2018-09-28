## install pacman and cmapR packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
if(!require(cmapR)){p_install('rhdf5');devtools::install_github("cmap/cmapR")}


# ###################################################################
# preprocess GCT file for subsequent ssGSEA/PTM-SEA analysis, e.g.
# - create site-centric GCT
# - create PTMsigDB compatible site identifier
# - create gene-centric GCT for ssGSEA; not implemented yet
preprocessGCT <- function(
  gct.str='',                         ## path to GCT file 
  level=c('site', 'gene'),            ## create  site or gene-centric reports; ('gene' not supported as of now)
  id.type=c('sm', 'wg'),              ## notation of site-ids: sm-Spectrum Mill; wg-Web Gestalt
  mode=c('mean', 'median', 'sd'),     ## how should multiple sites per gene be combine; sd - most variable (standard deviation) across sample columns 
  mod.res=c('S|T|Y', 'K'),            ## modified residue(s) 
  mod.type=c('-p', '-ac', '-ub'),     ## modification type
  acc.type=c('uniprot', 'refseq', 'symbol'),  ## 
  org=c('human', 'mouse', 'rat'),     ## supported organism; parameter not used right now; function 'RefSeq2UniProt' tries to determine organism automatically
  appenddim=T,                        ## see cmapR::write.gct()
  do.it=T                             ## flag, if FALSE nothing will be done; probably needed for to make this step optional in a FireCLoud WDL.
){
  
  ## imediatly return
  if(!do.it){
    return(0)
  }
  
  require(pacman)
  p_load(cmapR)
  p_load(magrittr)
  p_load(glue)
  
  #parse paramters----
  level <- match.arg(level)
  mode <- match.arg(mode)
  mod.res <- match.arg(mod.res)
  mod.type <- match.arg( mod.type )
  acc.type <- match.arg( acc.type )
  org <- match.arg(org)
  id.type <- match.arg(id.type)
  
  
  #import GCT---------
  gct <- try(parse.gctx(gct.str))
  
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
  
  #parse GCT object
  mat <- gct@mat
  n <- ncol(mat) ## number of data colums
  rid <- gct@rid
  if(is.null(rid)) rid <- rownames(mat)
  cid <- gct@cid
  if(is.null(cid)) cid <- colnames(mat)
  cdesc <- gct@cdesc
  #cdesc.ids <- colnames(cdesc)
  rdesc <- gct@rdesc
  #rdesc.ids <- colnames(rdesc)
  
  ## ################################
  ## create site-centric table
  if(level == 'site'){
    
    #row ids are site ids
    site.ids <- rid
    names(site.ids) <- rid
    site.ids <- gsub(' ', '', site.ids )
    
    # In case of RefSeq accession, remove all non-RefSeq accession numbers,
    # e.g. UniProt accessions of lab contaminants
    # then map RefSeq to UniProt
    if(acc.type == 'refseq'){
      #refseq only
      idx <- grep('^(NP_|XP_|YP_|ZP_)', site.ids)
      site.ids <- site.ids[idx]
      
      # update GCT
      rid <- rid[idx]
      if(n == 1){
        mat <- matrix(mat[idx, ], ncol=n)
        dimnames(mat) <- list(rid, cid)
      } else {
        mat <- mat[idx, ]
      }
      
      rdesc <- rdesc[idx, ]
      #rdesc <- matrix( rdesc[idx, ], ncol=length(rdesc.ids))
      #dimnames(rdesc) <- list(rid, rdesc.ids)
      
      #map to UniProt
      up <- RefSeq2UniProt( sub('^(.*?_.*?)_.*', '\\1', site.ids) )$id.map$id.mapped
      mapped.idx <- which(!is.na(up))
      
      ## update
      up <- up[mapped.idx]
      site.ids <- site.ids[mapped.idx]
      rid <- rid[mapped.idx]
      
      if(n == 1){
        mat <- matrix(mat[mapped.idx, ], ncol=ncol(mat))
        dimnames(mat) <- list(rid, cid)
      } else{
        mat <- mat[mapped.idx, ]
      }
      
      rdesc <- rdesc[mapped.idx, ]
      #rdesc <- matrix(rdesc[mapped.idx, ], ncol=length(rdesc.ids))
      #dimnames(rdesc) <- list(rid, rdesc.ids)
      
      site.ids <- paste(up, sub('^(.*?_.*?)(_.*)', '\\2', site.ids), sep='')
      names(site.ids) <- rid
      # use as ids
      #rownames(mat) <- rownames(rdesc) <- site.ids
      #rid <- site.ids
    }
    
    ## ##############################
    ## gene symbol + site
    ## EIF2S2.S2
    if(acc.type == 'symbol'){
      
      #map to UniProt
      up <- RefSeq2UniProt( sub('^(.*?)\\..*', '\\1', site.ids), keytype = 'SYMBOL' )$id.map$id.mapped
      mapped.idx <- which(!is.na(up))
      
      ## update
      up <- up[mapped.idx]
      site.ids <- site.ids[mapped.idx]
      rid <- rid[mapped.idx]
      
      if(n == 1){
        mat <- matrix(mat[mapped.idx, ], ncol=ncol(mat))
        dimnames(mat) <- list(rid, cid)
      } else {
        mat[mapped.idx, ]
      }
      
      rdesc <- rdesc[mapped.idx, ]
      #rdesc <- matrix(rdesc[mapped.idx, ], ncol=ncol(rdesc))
      #dimnames(rdesc) <- list(rid, rdesc.ids)
      
      ## creat PTM-SEA compatible ids
      site.ids <- glue("{up};{sub('^(.*?)\\\\.(.*)', '\\\\2', site.ids)}{mod.type}")
      names(site.ids) <- rid
    }
    ## ###############################
    ## Spectrum Mill VM site ids
    if(id.type == 'sm'){
      
      #localized sites----
      # - index of fully localized sites 
      loc.idx <- which( sapply (strsplit(sub('^.*_([1-9]_[0-9])_.*', '\\1', site.ids), '_'), function(x)length(unique(x)) ) == 1)
      
      # update GCT
      rid <- rid[ loc.idx ]
      site.ids <- site.ids[ loc.idx ]
      
      #mat <- mat[loc.idx, ]
      if(n == 1){
        mat <- matrix(mat[loc.idx, ], ncol=ncol(mat))
        dimnames(mat) <- list(rid, cid)
      } else {
        mat <- mat[loc.idx, ]
      }
      
      rdesc <- rdesc[loc.idx,]
      #rdesc <- matrix(rdesc[loc.idx, ], ncol=ncol(rdesc))
      #dimnames(rdesc) <- list(rid, rdesc.ids)
      
      ## variable sites
      site.var <- sub('^(.*?_.*?)_.*', '\\1', site.ids)
      names(site.var) <- names(site.ids)
      
      ## accession number plus modified residue
      all.sites <- lapply( strsplit(site.var, tolower(mod.res)), function(x){ 
        prot.acc=sub('_.*', '', x[1])
        x=sub('.*?_', '', x)
        paste(prot.acc, grep( toupper(mod.res), x, value=T), sep=';' )
        
      })
      
      #redundant sites on multiply modified peptides---- 
      ## Exclude sites on multiply phosphorylated peptides for which we
      ## have also detected singly phosphorylated version
      n.sites <- sapply(all.sites, length)
      all.sites.mult <- all.sites[which(n.sites > 1)]
      all.sites.sing <- all.sites[which(n.sites == 1)]
      all.sites.mult <- lapply(all.sites.mult, function(x){
        rm.idx=which(x %in% unlist(all.sites.sing))
        if(length(rm.idx)>0)
          x=x[-rm.idx]
        x
      })
      all.sites.mult <- all.sites.mult[ sapply(all.sites.mult, length) > 0 ]
      ## all sites as list
      all.sites <- append(all.sites.mult, all.sites.sing)
      all.sites.unique <- unique(unlist(all.sites))
      
      sites.dup <- unlist(all.sites)[ duplicated( unlist(all.sites)) ]
      sites.nondup <- setdiff(all.sites.unique, sites.dup)
      
      #prepare site-centric table----
      mat.ss <- matrix(NA, nrow=length(all.sites.unique), ncol=ncol(mat), dimnames=list(all.sites.unique, cid))
      
      ## TODO: parallelize loop
      # map of row indices between site-centric and original table
      #cl <- makeCluster(detectCores()-1)
      #registerDoParallel(cl)
      #map.idx <-  foreach(x = rownames(mat.ss)) %dopar% {function(x) grep(paste(x, '($|\\")', sep=''), all.sites)}
      #on.exit(stopCluster(cl))
      map.idx <- lapply(rownames(mat.ss), function(x) names(all.sites)[grep(paste(x, '($|\\")', sep=''), all.sites)])
      names(map.idx) <- rownames(mat.ss)
      
      #combine expression of duplicated sites, e.g. sites detected on multiple versions of multiply phosphorylated peptides
      for(s in sites.dup){
        if(n == 1){
          mat.tmp <- matrix(mat[map.idx[[s]], ], ncol=n)
        } else{
          mat.tmp <- mat[map.idx[[s]], ]
        }
        mat.ss[s, ] <- apply(mat.tmp, 2, eval(parse(text=mode)), na.rm=T)
      }
      mat.ss[names(map.idx[sites.nondup]), ] <- data.matrix( mat[unlist(map.idx[sites.nondup]), ] )
      
      #update GCT----
      ## - multiple sites (rids) have been combined.
      ## - pick single rid for combined sites
      idx <- sapply(map.idx, function(x)x[1])
      #rid.org <- rid
      #rdesc.org <- rdesc
      rid <- rid[ idx ]
      names(rid) <- names(map.idx)
      rdesc <- rdesc[idx, ]
      VMsitesAll <- sapply(map.idx, function(x) paste(names(all.sites)[x], collapse='|'))
      rdesc <- data.frame(rdesc, VMsitesAll)
      
      #create site ids-----
      ptm.site.ids <- paste(names(map.idx), mod.type, sep=ifelse(length(grep('^-', mod.type)) > 0, '','-'))
      
    } else { # end if Spectrum Mill VM ids
      ## assume input is site centric
      ptm.site.ids <- site.ids
      mat.ss <- mat
    }
    
    ##rid <- ptm.site.ids
    rownames(mat.ss) <- rownames(rdesc) <- rid <- ptm.site.ids
    #export GCT
    gct@mat <- mat.ss
    gct@rid <- rid
    gct@cid <- cid
    gct@cdesc <- cdesc
    gct@rdesc <- rdesc
    
    fn <- ifelse(appenddim,
                 paste(sub('_n[0-9]*x[0-9]*\\.gct$', paste('_',level, '-centric', sep=''), gct.str), '', sep=''),
                 paste(sub('\\.gct',paste('_',level, '-centric.gct', sep=''), gct.str), '', sep='')
    ) 
    write.gct(gct, fn, appenddim = appenddim)
    
  }
  
  return(0)
}

## ########################################################################
## create UniProt-centric accesion numbers
## -map RefSeq accessions or gene symbol to UniProt accession
RefSeq2UniProt <- function(ids,                          ## character vector of accessions
                           n.try=10,                     ## maximal number of accessions to test in order to dermine the source organism
                           keytype=c('REFSEQ', 'SYMBOL') ## RefSeq or gene symbols in 'ids'  
                           ){
  
  require(pacman)
  p_load(RSQLite)
  p_load(org.Hs.eg.db)
  p_load(org.Mm.eg.db)
  p_load(org.Rn.eg.db)
  p_load(org.Dr.eg.db)
  
  ## ###################################
  ##           id type
  ## ###################################
  keytype <- match.arg(keytype)
  
  ## ###################################
  ##        extract query strings
  ## ###################################
  id.query <- sub('(\\.|;).*', '', ids) ## first id
  names(id.query) <- ids
  
  ## ###################################
  ##          determine organism
  ## ###################################
  orgtype <- 'UNKNOWN'
  
  # try human 
  if(orgtype == 'UNKNOWN'){ 
    id.map.tmp <- try( mapIds(org.Hs.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('UNIPROT'), keytype=keytype, multiVals='first') )
    if(class(id.map.tmp) != 'try-error'){
      orgtype='HSA'
    }
  }
  # try mouse
  if(orgtype == 'UNKNOWN'){ 
    id.map.tmp <- try( mapIds(org.Mm.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('UNIPROT'), keytype=keytype, multiVals='first') )
    if(class(id.map.tmp) != 'try-error'){
      orgtype='MMU'
    }
  }
  # try rat
  if(orgtype == 'UNKNOWN'){ 
    id.map.tmp <- try( mapIds(org.Rn.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('UNIPROT'), keytype=keytype, multiVals='first') )
    if(class(id.map.tmp) != 'try-error'){
      orgtype='RNO'
    }
  }
  # try zebrafish
  if(orgtype == 'UNKNOWN'){ 
    id.map.tmp <- try( mapIds(org.Dr.eg.db, keys=id.query[ sample(1:length(id.query), n.try)] , column=c('UNIPROT'), keytype=keytype, multiVals='first') )
    if(class(id.map.tmp) != 'try-error'){
      orgtype='DRE'
    }
  }
  
  ## ##################################
  ## map
  if( orgtype != 'UNKNOWN'){
    if(orgtype == 'HSA')
      id.map.tmp <- try(mapIds(org.Hs.eg.db, keys=id.query , column=c('UNIPROT'), keytype=keytype, multiVals='first'))
    if(orgtype == 'MMU')
      id.map.tmp <- try(mapIds(org.Mm.eg.db, keys=id.query , column=c('UNIPROT'), keytype=keytype, multiVals='first'))
    if(orgtype == 'RNO')
      id.map.tmp <- try(mapIds(org.Rn.eg.db, keys=id.query , column=c('UNIPROT'), keytype=keytype, multiVals='first'))
    if(orgtype == 'DRE')
      id.map.tmp <- try(mapIds(org.Dr.eg.db, keys=id.query , column=c('UNIPROT'), keytype=keytype, multiVals='first'))
  } else {
    id.map.tmp <- c()
  }
  
  if(class(id.map.tmp) == 'try-error' | is.null( class(id.map.tmp) ) | class(id.map.tmp) == 'NULL' ){
    
    id.map <- data.frame(id=names(id.query), id.query=id.query, id.mapped=names(id.query), id.concat=ids, stringsAsFactors=F)

    keytype <- 'UNKNOWN'
    
  } else {
    
    #id.map.tmp[which(is.na(id.map.tmp))] <- 'NotFound'
    id.map <- data.frame(id=names(id.query), id.query=id.query, id.mapped=id.map.tmp, id.concat=paste(ids, id.map.tmp, sep='_'), stringsAsFactors=F)
  }
  
  ## results
  res <- list()
  res[[1]] <- keytype
  res[[2]] <- id.map
  res[[3]] <- orgtype
  names(res) <- c('keytype', 'id.map', 'orgtype')
  
  return(res)
}

## ##########################################################################
## function to convert an output file from PAM-based marker selection
## into a GCT 1.3 file as input for PTM-SEA
prepare_pam_output <- function( csv.str,  ## path to input csv
                                expr='^(contrast-|Fold Change)',
                                ofile=sub('\\.csv', '\\.gct', sub('.*/','', csv.str))
){
  require(pacman)
  p_load(glue)
  
  ## import
  csv <- read_csv(csv.str) %>% as.data.frame
  
  expr.idx <- grep(expr, colnames(csv))
  
  if(length(expr.idx) == 0) stop(glue('\n\nNo columns matching "{expr}" found in header of input file {csv.str}.\n\n'))
  
  
  mat <- csv[, expr.idx]  
  cid <- colnames(csv)[expr.idx]
  rdesc <- csv[, setdiff( colnames(csv), cid)]
  rid <- csv[, 1] 
  
  ## assemble GCT
  gct <- new('GCT')
  gct@mat <- data.matrix(mat)
  gct@rdesc <- rdesc
  gct@rid <- rid
  gct@cid <- cid
  
  fn <- sub('.*/', '', csv.str) %>% sub('\\.csv', '\\.gct', .)
  
  write.gct(gct, ofile = ofile, appenddim = F)
  
  return(0)
}


