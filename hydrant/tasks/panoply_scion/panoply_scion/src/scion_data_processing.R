#pre-processes data files for SC-ION
#allows for random permutations if desired

scion_data_processing <- function(pome.gct.file, mrna.gct.file, TF, prefix,permute=NULL,dim="row",na.max=0.7,TF.filter=T,type=c("SM","FP")){
  #need to change the seed for each permutation. otherwise, the shuffling is exactly the same when running in parallel
  if(!is.null(permute)){
    set.seed(permute)
  }
  
  #read in matrices
  reg <- parse_gctx(pome.gct.file)
  reg.mat <- reg@mat
  target <- parse_gctx(mrna.gct.file)
  target.mat <- target@mat
  
  #make sure the columns are in the same order in both matrices
  reg.mat <- reg.mat[,order(colnames(reg.mat))]
  target.mat <- target.mat[,order(colnames(target.mat))]
  #and, make sure the columns are the same
  reg.mat <- reg.mat[,colnames(reg.mat)%in%colnames(target.mat)]
  target.mat <- target.mat[,colnames(target.mat)%in%colnames(reg.mat)]
  
  #target data do not need to be normalized - they are normalized before clustering/network inference by default
  #target.mat <- rowNorms(target.mat, type="z")
  
  #permute data matrix if desired
  if(!is.null(permute)){
    #apply missing value filter BEFORE permuation!
    target.mat <- target.mat[(rowSums(is.na(target.mat))/dim(target.mat)[2])<=na.max,]
    if(dim=="col"){
      #permute the target matrix column-wise (by sample)
      target.names <- row.names(target.mat)
      target.mat <- apply(target.mat,2,sample)
      row.names(target.mat) <- target.names
    }else{
      #permute the target matrix row-wise (by gene)
      target.mat <- t(apply(target.mat,1,sample))
    }
  }
  
  #we need to append group ID or site information to the gene symbol
  if(prefix %in% c("phosphoproteome","acetylome","ubiquitylome")){
    #PTMs
    sites <- lapply(row.names(reg.mat),function(x){
      mysites = unlist(strsplit(x,"_"))
      #for SpectrumMill, the sites always start with an uppercase letter and end with a lowercase letter. this seems to be unique to the sites and works for now but may need revisited.
      if(type=="SM"){
        mysite <- mysites[grepl("[[:lower:]]$",mysites) & grepl("^[[:upper:]]",mysites)]
      }
      #FragPipe has site ID always at the end
      if(type=="FP"){
        mysite <- mysites[length(mysites)]
      }
      mysite <- str_trim(mysite)
      return(mysite)
    })
    sites <- as.character(sites)
    new.symbols <- paste(reg@rdesc$geneSymbol,sites,sep="_")
    row.names(reg.mat)=new.symbols
    
    #remove rows which don't have a mapped site
    #this sometimes happen if we try and convert sites for other datasets
    #for SM only
    if(type=="SM"){
      reg.mat <- reg.mat[grepl('[[:lower:]]',row.names(reg.mat)) & !grepl('character',row.names(reg.mat)),]
    }
    
    #make a list of the TFs with site information
    #only need to do this if filtering for TFs
    if(TF.filter){
      TF.list <- NULL
      for(i in 1:dim(TF)[1]){
        if(TF[i,1]%in%reg@rdesc$geneSymbol){
          TF.list = c(TF.list,row.names(reg.mat)[reg@rdesc$geneSymbol%in%TF[i,1]])
        }
      }
      TF.list <- unique(TF.list)
    }
    
    #permute data matrix if desired
    if(!is.null(permute)){
      #filter missing values BEFORE permuting!
      reg.mat <- reg.mat[(rowSums(is.na(reg.mat))/dim(reg.mat)[2])<=na.max,]
      if(dim=="col"){
        #permute the regulator matrix col-wise (by sample)
        reg.names <- row.names(reg.mat)
        reg.mat <- apply(reg.mat,2,sample)
        row.names(reg.mat) <- reg.names
      }else{
        #permute the regulator matrix row-wise (by site)
        reg.mat <- t(apply(reg.mat,1,sample))
      }
    }
    
    #filter the regulator matrix for the TFs (if desired)
    if(TF.filter){
      reg.mat <- reg.mat[row.names(reg.mat)%in%TF.list,]
    }
    #ensure there are no duplicated row names
    reg.mat <- reg.mat[!duplicated(row.names(reg.mat)),]
    
    #make the clustering matrix
    #for PTMs we just stick the two together
    cluster.mat = rbind(reg.mat,target.mat)
    cluster.mat = unique(cluster.mat)
  }else{
    #aggregate expression for duplicated gene symbols
    if(sum(duplicated(reg@rdesc$geneSymbol)) > 0){
      mat <- aggregate(reg.mat, list(reg@rdesc$geneSymbol), function(x) mean(x, na.rm=T))
      #remove missing gene symbols
      mat <- mat[mat$Group.1!="",]
      rid <- mat$Group.1
      mat <- mat[, -c(1)]
      row.names(mat) <- rid
    }
    reg.mat <- mat
    
    #permute data matrix if desired
    if(!is.null(permute)){
      #filter missing values BEFORE permuting!
      reg.mat <- reg.mat[(rowSums(is.na(reg.mat))/dim(reg.mat)[2])<=na.max,]
      if(dim=="col"){
        #permute the regulator matrix col-wise (by sample)
        reg.names <- row.names(reg.mat)
        reg.mat <- apply(reg.mat,2,sample)
        row.names(reg.mat) <- reg.names
      }else{
        #permute the regulator matrix row-wise (by protein)
        reg.mat <- t(apply(reg.mat,1,sample))
      }
    }
    
    #filter the regulator matrix for the TFs (if desired)
    if(TF.filter){
      TF.list <- TF[,1]
      reg.mat <- reg.mat[row.names(reg.mat)%in%TF.list,]
    }
    
    #make the clustering matrix
    #in the clustering matrix, we use regulator data for regulators and target data for targets
    target.mat.noTFs = target.mat[!row.names(target.mat)%in%row.names(reg.mat),]
    cluster.mat = rbind(reg.mat,target.mat.noTFs)
  }
  
  #missing value filter (default 70%)
  #for non-permutation this is the only filter
  #for permutation this filters any features that exceed na.max after permutation
  target.mat <- target.mat[(rowSums(is.na(target.mat))/dim(target.mat)[2])<=na.max,]
  reg.mat <- reg.mat[(rowSums(is.na(reg.mat))/dim(reg.mat)[2])<=na.max,]
  cluster.mat <- cluster.mat[(rowSums(is.na(cluster.mat))/dim(cluster.mat)[2])<=na.max,]
  
  #replace NA with zeros
  target.mat[is.na(target.mat)]<-0
  reg.mat[is.na(reg.mat)]<-0
  cluster.mat[is.na(cluster.mat)]<-0
  
  #return the regulator, target, and clustering matrices
  return(list("reg"=reg.mat,"target"=target.mat,"cluster"=cluster.mat))
}