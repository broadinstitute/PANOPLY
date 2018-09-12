






## ##########################################################################
## function to convert an output file from PAM-based marker selection
## into a GCT 1.3 file as input for PTM-SEA
prepare_pam_output <- function( csv.str,  ## path to input csv
                                expr='^contrast-',
                                ofile=sub('\\.csv', '\\.gct', sub('.*/','', csv.str))
){
  
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
}



## ########################################################################
## function to map RefSeq accession numbers  to UniProt accession numbers
RefSeq2UniProt <- function(ids, n.try=10 ){
  
  ## id mapping
  p_load(RSQLite)
  p_load(org.Hs.eg.db)
  p_load(org.Mm.eg.db)
  p_load(org.Rn.eg.db)
  p_load(org.Dr.eg.db)
  
  ## ###################################
  ##           id type
  ## ###################################
  keytype <- 'REFSEQ'
  
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
  
  #cat(orgtype)
  
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
    #cat('test2\n')
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


# ###################################################################
#preprocess GCT file for subsequent ssGSEA/ssPSEA analysis, e.g.
# - create site-centric GCT
# - map RefSeq to UniProt
# - create gene-centric GCT
preprocessGCT <- function(#gct.str='/media/sf_Karsten/Projects/20180122_CPATC3_PTRC/Westbrook/data/phosphoproteome-ratio-norm-NArm.gct',
  gct.str='//flynn-cifs/prot_proteomics/LabMembers/Karsten/Projects/20180122_CPATC3_PTRC/Westbrook/data/phosphoproteome-ratio-norm-NArm.gct', 
  level=c('site', 'gene'), 
  mode=c('mean', 'median', 'sd'),
  mod.res=c('S|T|Y', 'K'), #modified residue(s)
  mod.type=c('-p', '-ac', '-ub'),
  acc.type=c('uniprot', 'refseq'),
  org=c('human', 'mouse', 'rat'),
  appenddim=T,
  do.it=T   ## flag, if FALSE nothing will be done 
){
  
  ## imediatly return
  if(!do.it){
    return(0)
  }
  
  require(pacman)
  p_load(cmapR)
  p_load(magrittr)
  
  #parse paramters----
  level <- match.arg(level)
  mode <- match.arg(mode)
  mod.res <- match.arg(mod.res)
  mod.type <- match.arg( mod.type )
  acc.type <- match.arg( acc.type )
  org <- match.arg(org)
  
  
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
  rid <- gct@rid
  cid <- gct@cid
  cdesc <- gct@cdesc
  rdesc <- gct@rdesc
  
  #create site-centric table-------------------------------------------------
  ## assuming Spectrum Mill Site ids, e.g: NP_001611.1_S4386s _1_1_4386_4386
  if(level == 'site'){
    
    #row ids are site ids
    site.ids <- rid
    names(site.ids) <- rid
    site.ids <- gsub(' ', '', site.ids )
    
    # In case of RefSeq database, remove all non-RefSeq accession numbers,
    # e.g. UniProt accessions of lab contaminants
    # then map RefSeq to UniProt
    if(acc.type == 'refseq'){
      #refseq only
      idx <- grep('^(NP_|XP_|YP_|ZP_)', site.ids)
      site.ids <- site.ids[idx]
      
      # update GCT
      mat <- mat[idx, ]
      rid <- rid[idx]
      rdesc <- rdesc[idx, ]
      
      #map to UniProt
      up <- RefSeq2UniProt( sub('^(.*?_.*?)_.*', '\\1', site.ids) )$id.map$id.mapped
      mapped.idx <- which(!is.na(up))
      
      ## update
      up <- up[mapped.idx]
      site.ids <- site.ids[mapped.idx]
      rid <- rid[mapped.idx]
      mat <- mat[mapped.idx, ]
      rdesc <- rdesc[mapped.idx, ]
      
      site.ids <- paste(up, sub('^(.*?_.*?)(_.*)', '\\2', site.ids), sep='')
      names(site.ids) <- rid
      # use as ids
      #rownames(mat) <- rownames(rdesc) <- site.ids
      #rid <- site.ids
    }
    
    #localized sites----
    # - index of fully localized sites 
    loc.idx <- which( sapply (strsplit(sub('^.*_([1-9]_[0-9])_.*', '\\1', site.ids), '_'), function(x)length(unique(x)) ) == 1)
    
    # update GCT
    rid <- rid[ loc.idx ]
    site.ids <- site.ids[ loc.idx ]
    mat <- mat[loc.idx, ]
    rdesc <- rdesc[loc.idx,]
    
    #parse UniProt accession----
    ## the code below works for UniProt accession numbers. Leave as is for the time being
    ## and merge at a later point
    #if(acc.type == 'uniprot'){
    
    ## variable sites
    site.var <- sub('^(.*?_.*?)_.*', '\\1', site.ids)
    names(site.var) <- names(site.ids)
    
    ## accession number plus modified residue
    all.sites <- lapply( strsplit(site.var, tolower(mod.res)), function(x){ 
      prot.acc=sub('_.*', '', x[1])
      x=sub('.*?_', '', x)
      paste(prot.acc, grep( toupper(mod.res), x, value=T), sep=';' )
      
    })
    #names(all.sites) <- site.ids
    #}
    #parse RefSeq accessions----
    ## RefSeq accesion numbers
    # if(acc.type == 'refseq'){
    #   
    #   ## variable sites
    #   site.var <- sub('^(.*?_.*?_.*?)_.*', '\\1', site.ids)
    #   
    #   ## accession number plus modified residue
    #   all.sites <- lapply( strsplit(site.var, '(s|t|y)'), function(x){ 
    #     prot.acc=sub('_[S|T|Y].*', '', x[1])
    #     x=sub('^.*_([S|T|Y].*)$', '\\1', x)
    #     paste(prot.acc, grep( '(S|T|Y)', x, value=T), sep=';' )
    #     
    # #  })
    #  names(all.sites) <- site.ids
    #}
    
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
    
    # map of row indices between site-centric and original table
    ##map.idx <- lapply(rownames(mat.ss), function(x) grep(x, all.sites))
    #cat('yeppp111\n')
    #cl <- makeCluster(detectCores()-1)
    #registerDoParallel(cl)
    #map.idx <-  foreach(x = rownames(mat.ss)) %dopar% {function(x) grep(paste(x, '($|\\")', sep=''), all.sites)}
    #on.exit(stopCluster(cl))
    #cat('yeppp\n')
    map.idx <- lapply(rownames(mat.ss), function(x) names(all.sites)[grep(paste(x, '($|\\")', sep=''), all.sites)])
    names(map.idx) <- rownames(mat.ss)
    
    #combine expression of duplicated sites
    for(s in sites.dup){
      mat.tmp <- mat[map.idx[[s]], ]
      mat.ss[s, ] <- apply(mat.tmp, 2, eval(parse(text=mode)), na.rm=T)
    }
    mat.ss[names(map.idx[sites.nondup]), ] <- data.matrix( mat[unlist(map.idx[sites.nondup]), ] )
    
    
    #mat[names(all.sites[map.idx[[sites.dup[1]]]]), ]
    
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
    
    return(0)
  }
  
}

