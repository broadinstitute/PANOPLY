#' SCION implementation to run on Terra. Combines preprocessing and network inference into one function
#' @author Natalie M Clark
#' @param prefix ome type of the regulator (proteome, phosphoproteome, acetylome, ubiquitylome)
#' @param pome.gct.file path to the regulator (proteome) GCT file
#' @param mrna.gct.file path to the target (RNA) GCT file
#' @param TF.file path to the file containing a list of gene symbols that correspond to known transcription factors (TFs). It is assumed this file contains one gene symbol per line and has no header.
#' @param permute number of random permutations to perform for edge trimming. Default NULL, which means no permutations will be performed
#' @param dim dimension on which to permute. options are "row" or "col". default "col"
#' @param na.max minimum proportion of missing values between 0 and 1. all features with proportion of missing values < na.max are filtered. any remaining missing values will be imputed with 0. default 0 (use only fully quantified features). 
#' @param cluster boolean option to cluster using k-means prior to network inference. default TRUE
#' @param type type of proteomics data, i.e. the engine used to search the proteomics data. options are "SM" (SpectrumMill) or "FP" (Fragpipe). default "SM"
#' @param dir.name name of directory to save results. default "exp". if using permutations, this parameter is ignored, and the directories are named after the permutation number (1,2,3...)
#' @param weightthreshold threshold for edge trimming. all edges with weight < weightthreshold are removed from the network. minimum value of 0. default 0.
#' @param normalize boolean to normalize edge weights to a [0,1] scale. default FALSE
#' @param num.cores number of cores to use for parallelization. num.cores-1 will be used for parallelization. if num.cores < 3, no parallelization is performed. default 1.
#' @param connect.hubs boolean to connect the hubs between clusters. this parameter is ignored if clustering is not performed. default TRUE
#' @param verbose boolean to display detailed output. default FALSE
SCION <- function (prefix, pome.gct.file, mrna.gct.file, TF.file, permute=NULL, dim="col", na.max=0, cluster=T,dir.name="exp",type="SM", weightthreshold=0, normalize=FALSE,num.cores=1,connect.hubs=T,verbose=F) {
  
  # Source the files that we need
  source("scion_data_processing.R") #pre-processing data tables
  source("kmeans_clustering.R") #clustering
  source("RS.Get.Weight.Matrix.parallel.R") #network inference
  
  #load packages
  library(cmapR)
  library(stringr)
  library(dplyr)
  library(cluster)
  library(randomForest)
  library(doParallel)
  
  #create directory to save results
  if(!is.null(permute)){
    my.dir <- as.character(permute)
    cat(paste("Starting permutation ",permute,"\n",sep=""))
  }else{
    my.dir <- dir.name
  }
  
  if(!dir.exists(my.dir)){
    dir.create(my.dir)
  }
  
  #first need to pre-process the data
  cat("Processing data tables\n")
  #read in the TF file. 1 TF per line, no header.
  TF <- read.table(TF.file, header=FALSE)
  tables <- scion_data_processing(pome.gct.file, mrna.gct.file, TF, prefix,permute,dim,na.max,type=type) 
  
  #assign pre-processed data tables
  mytargetdata <- as.data.frame(tables$target)
  myregdata <- as.data.frame(tables$reg)
  myclustdata <- as.data.frame(tables$cluster)
  
  #save data tables
  write.csv(mytargetdata,paste(my.dir,"/target_mat.csv",sep=""))
  write.csv(myregdata,paste(my.dir,"/reg_mat.csv",sep=""))
  write.csv(myclustdata,paste(my.dir,"/clust_mat.csv",sep=""))
  
  #cluster using kmeans (if desired)
  if(cluster==T){
    cat("Clustering\n")
    #choose starting value of k
    #if there are many more regulators than samples, scale starting k based on regulators
    if(dim(myregdata)[1]>10*dim(myregdata)[2]){
      kmid <- dim(myregdata)[1]*0.1
    }else{
      #otherwise scale k based on the smallest dimension (regulators or samples)
      kmid <- min(dim(myregdata)[1],dim(myregdata)[2])
    }
    clusterresults <- kmeans_clustering(myclustdata, kmid)
    
    #save the results
    write.csv(clusterresults,paste(my.dir,"/clusters.csv",sep=""))
  }
  
  #infer a network using GENIE3 on each cluster
  SCION_infer <- function(mytargetdata,myregdata,clusterresults,weightthreshold,prefix,permute=NULL,normalize=FALSE,verbose=F){
    
    finalnetwork <- data.frame(Regulator=character(), Interaction=character(), Target=character(), Weight=double(), Cluster=numeric(),
                               stringsAsFactors=FALSE)
    myhubs=NULL
    cat("Inferring cluster-specific networks.\n")
    for (i in 1:max(clusterresults$clusters)){
      mygenes = row.names(clusterresults)[clusterresults$clusters==i]
      clustertargetdata = mytargetdata[row.names(mytargetdata)%in%mygenes,]
      clusterregdata = myregdata[row.names(myregdata)%in%mygenes,]
      
      #we need at least one target and at least one regulator
      if (dim(clustertargetdata)[1]<1 | dim(clusterregdata)[1]<1){
        next
      }
      
      #infer the network
      if(verbose){
        cat(paste("Inferring network for cluster ",i,
                  " with ", dim(clusterregdata)[1], " regulators and ",
                  dim(clustertargetdata)[1], " targets.\n",sep=""))
      }
      network = RS.Get.Weight.Matrix(t(clustertargetdata),t(clusterregdata),normalize=normalize,num.cores=num.cores)
      #if network inference failed, move on
      if (is.null(network)){
        next
      }
      
      #now make a new network where we eliminate all the low confidence edges
      trimmednet = data.frame(ifelse(network<weightthreshold,NaN,network))
      
      #translate the trimmed network into a table we can import into cytoscape
      networktable = data.frame(Regulator=character(), Interaction=character(), Target=character(), Weight=double(),Cluster=numeric(),
                                stringsAsFactors=FALSE)
      row=1
      for (j in 1:dim(trimmednet)[1]){
        for (k in 1:dim(trimmednet)[2]){
          #skip NaNs as these have no edge
          if (is.na(trimmednet[j,k])){
            next
          }else{
            networktable[row,] = cbind(colnames(trimmednet)[k],"regulates",rownames(trimmednet)[j],trimmednet[j,k],i)
            row = row+1
          }
        }
      }
      
      #save this network
      finalnetwork = rbind(finalnetwork,networktable)
      
      #get the hub gene (highest outdegree) and save it to connect the clusters later
      #if there is a tie, we save both hubs
      myregs = unique(networktable$Regulator)
      for (j in 1:length(myregs)){
        numedges = sum(networktable$Regulator%in%myregs[j])
        if (j==1){
          hubedges = numedges
          hub = myregs[j]
        }else if (numedges>hubedges){
          hubedges = numedges
          hub = myregs[j]
        }else if (numedges==hubedges){
          hub = rbind(hub,myregs[j])
        }
      }
      if (!exists("myhubs")){
        myhubs = hub
      } else{
        myhubs = rbind(myhubs,hub)
      }
    }
    
    #connect the hubs for each cluster
    if (connect.hubs && exists("myhubs") && length(myhubs)>2){
      cat("Connecting hubs\n")
      hubtargetdata = mytargetdata[row.names(mytargetdata)%in%myhubs,]
      hubregdata = myregdata[row.names(myregdata)%in%myhubs,]
      if (dim(hubtargetdata)[1]==0){
        #strip the PTM information so that we can get the targets
        genes <- unlist(strsplit(myhubs,'_'))
        genes <- genes[seq(1,length(genes),by=2)]
        hubtargetdata = mytargetdata[row.names(mytargetdata)%in%genes,]
      }
      network = RS.Get.Weight.Matrix(t(hubtargetdata),t(hubregdata),normalize=normalize)
      
      #now make a new network where we eliminate all the low confidence edges
      trimmednet = data.frame(ifelse(network<weightthreshold,NaN,network))
      
      #translate the trimmed network into a table we can import into cytoscape
      networktable = data.frame(Regulator=character(), Interaction=character(), Target=character(), Weight=double(),Cluster=numeric(),
                                stringsAsFactors=FALSE)
      row=1
      for (j in 1:dim(trimmednet)[1]){
        for (k in 1:dim(trimmednet)[2]){
          #skip NaNs as these have no edge
          if (is.na(trimmednet[j,k])){
            next
          }else{
            networktable[row,] = cbind(colnames(trimmednet)[k],"regulates",rownames(trimmednet)[j],trimmednet[j,k],"")
            row = row+1
          }
        }
      }
      #save this network
      finalnetwork = rbind(finalnetwork,networktable)
    }
    
    
    #for the PTMs, strip the site information and save in a separate column
    if(prefix %in% c('acetylome','phosphoproteome','ubiquitylome')){
      sites <- unlist(strsplit(finalnetwork$Regulator,"_"))
      #this doesn't always work - if it doesn't, just skip the rest
      if(length(sites[seq(1,length(sites)-1,by=2)]) == length(finalnetwork$Regulator)){
        finalnetwork$Regulator <- sites[seq(1,length(sites)-1,by=2)]
        finalnetwork <- add_column(finalnetwork, sites[seq(2,length(sites),by=2)],.after="Regulator")
        colnames(finalnetwork)[2]="Site"
      }
    }
    
    #save the final network and the clustering information
    if(!is.null(permute)){
      cat(paste("Saving network for permutation ",permute,"\n",sep=""))
      write.table(finalnetwork,paste (my.dir,"/",prefix, '-SCION-network-permutation-',permute,'.tsv', sep=''),row.names=FALSE,quote=FALSE,sep='\t')
    }else{
      write.table(finalnetwork,paste (my.dir, "/", prefix, '-SCION-network-full.tsv', sep=''),row.names=FALSE,quote=FALSE,sep='\t')
    }
  }
  
  #infer the network
  #set seed so network inference is the same each time
  set.seed(2023)
  SCION_infer(mytargetdata,myregdata,clusterresults,weightthreshold,prefix,permute,normalize,verbose)
  cat("Done\n")
}