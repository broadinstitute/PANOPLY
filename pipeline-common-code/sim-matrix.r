
## install and load libraries automatically
# library (NMF)
if (!require("pacman")) install.packages ("pacman")
pacman::p_load (NMF)

nmf.options (grid.patch=TRUE)


# 
# jaccard.similarity <- function (S, T) {
#   # See Ioffe, S. (2010) Improved consistent sampling, weighted minhash and l1 sketching. 
#   #     10th IEEE Conference on Data Mining (ICDM), 246–255.
#   S1 <- norm (as.matrix (S), type='1')
#   T1 <- norm (as.matrix (T), type='1')
#   S_T1 <- norm (as.matrix (S-T), type='1')
#   
#   j.ST <- (S1 + T1 - S_T1) / (S1 + T1 + S_T1)
#   return (j.ST)
# }
# 
# 
# marker.similarity.matrix <- function (marker.lists, output.prefix) {
#   # calculate Jaccard similarity matrix 
#   # for all possible pairs of markers in the marker.lists list
#   # the names of the lists in marker.lists is used in the output
#   # each item in the marker.list must be a dataset with 2 columns:
#   # id (protein/gene name, etc.) and 0-1 measure (p-value, weight, etc.)
#   # if a marker is not present in a list when comparing a pair of lists,
#   # that marker is assigned a measure of 1
#   
#   nlists <- length (marker.lists)
#   # initialize similarity matrix
#   sim.matrix <- diag (nlists)
#   colnames (sim.matrix) <- rownames (sim.matrix) <- names (marker.lists)
#   
#   # calculate similarity
#   for (i in 2:nlists)
#     for (j in 1:(i-1)) {
#       markers <- merge (marker.lists[[i]], marker.lists[[j]], by=1, all=TRUE)
#       markers [is.na (markers[,2]), 2] <- 1
#       markers [is.na (markers[,3]), 3] <- 1
#       sim.matrix [i,j] <- jaccard.similarity (markers[,2], markers[,3])
#       sim.matrix [j,i] <- sim.matrix [i,j]
#     }
#   write.csv (sim.matrix, paste (output.prefix, '.csv', sep=''))
#   
#   # plot results
#   pdf (paste (output.prefix, '.pdf', sep=''), width=8, height=8)
#   aheatmap (sim.matrix, color='-heat')
#   dev.off ()
# }
# 
# 

jaccard.similarity <- function (S, T) {
  j.ST <- length (intersect (S, T)) / length (union (S, T))
  return (j.ST)
}


similarity.matrix <- function (marker.sets, output.prefix) {
  nlists <- length (marker.sets)
  # initialize similarity matrix
  sim.matrix <- diag (nlists)
  colnames (sim.matrix) <- rownames (sim.matrix) <- names (marker.sets)
  
  # calculate similarity
  for (i in 2:nlists)
    for (j in 1:(i-1)) {
      sim.matrix [i,j] <- jaccard.similarity (marker.sets[[i]], marker.sets[[j]])
      sim.matrix [j,i] <- sim.matrix [i,j]
    }
  write.csv (sim.matrix, paste (output.prefix, '.csv', sep=''))
  
  # plot results
  pdf (paste (output.prefix, '.pdf', sep=''), width=8, height=8)
  aheatmap (sim.matrix, color='YlOrRd:100', scale='none', border_color='grey90', 
            breaks=seq(0,1,0.01))
  dev.off ()
}



similarity <- function (file.list, output.prefix, annotation.list=NULL, marker.col=1, fc.col=6) {
  # calculate and plot similarity or markser selection output file in file.list
  # thr first column in the file is assumed to contain the list of markers
  # annotation list provides names to use for each file.list (defaults to file name)
  markers <- NULL
  annotation <- NULL
  for (i in 1:length(file.list)) {
    f <- file.list[i]
	annot <- ifelse (is.null (annotation.list), f, annotation.list[i])
    d <- read.csv (f)
	fc.up <- d[,fc.col] > 1
	fc.down <- d[,fc.col] < 1
	if (sum (fc.up) > 1) {
	  markers <- c (markers, list (d [fc.up, marker.col]))
	  annotation <- c (annotation, paste (annot, '.FC+', sep=''))
	}
	if (sum (fc.down) > 1) {
	  markers <- c (markers, list (d [fc.down, marker.col]))
	  annotation <- c (annotation, paste (annot, '.FC-', sep=''))
	}
  }
  names (markers) <- annotation
  
  similarity.matrix (markers, output.prefix)
}


##
## marker lists included in similarity matrix
## (collected from multiple sources/analyses/locations)
##
master.list <- NULL
master.annotation <- NULL

add.to.list <- function (files, annot) {
  present <- file.exists (files)
  master.list <<- c (master.list, files[present])
  master.annotation <<- c (master.annotation, annot[present])
}

                                       

## directories to compare
dir.list <- c ('current-tcga-proteome', 'current-whim-proteome', 'tcga+whim/combined-proteome-v1')
type.list <- c ('tcga', 'whim', 'tcga+whim')


for (i in 1:length(dir.list)) {
  d <- dir.list[i]
  type <- type.list[i]
  
  ## gene mutations
  mut.dir <- file.path ('..', d, 'mutation')
  genes <- c ('GATA3', 'PIK3CA', 'TP53')
  files <- file.path (mut.dir, paste ('gene-', genes, '/', genes, '-markers-fdr0.01.csv', sep=''))
  annotation <-  paste ('Mutation', genes, type, sep='-')
  add.to.list (files, annotation)


  ## clustering
  clus.dir <- file.path ('..', d, 'clustering')
  classes <- c (1, 2, 3)
  files <- file.path (clus.dir, paste ('cluster-class.', classes, '-markers-fdr0.01.csv', sep=''))
  annotation <- paste ('Cluster', classes, type, sep='-')
  add.to.list (files, annotation)
  # clustering when whim samples are included: only 2 clusters present
  files <- file.path (clus.dir, c ('whim-cluster-auto-analysis-markers-fdr0.01.csv',
								   'tcga+whim-cluster-auto-analysis-markers-fdr0.01.csv'))
  annotation <- c ('Cluster-1-whim', 'Cluster-1-tcga+whim')
  add.to.list (files, annotation)  
                        

  ## pam50
  pam50.dir <- file.path ('..', d, 'pam50')
  groups <- c ('Basal', 'Her2', 'LumA', 'LumB')
  files <- file.path (pam50.dir, paste ('pam50only-class.', groups, '-markers-fdr0.01.csv', sep=''))
  annotation <- paste ('PAM50', groups, type, sep='-')
  add.to.list (files, annotation)

  groups2 <- c ('lumA_v_lumB', 'basal_v_luminal')
  files <- file.path (pam50.dir, paste (groups2, '-markers-fdr0.01.csv', sep=''))
  annotation <- paste (groups2, type, sep='-')
  add.to.list (files, annotation)

  # ## association
  # assoc.dir <- file.path ('..', d, 'association')
  # groups <- c ('er', 'pr')
  # files <- file.path (assoc.dir, paste (groups, '-analysis-markers-fdr0.01.csv', sep=''))
  # annotation <- paste (groups, type, sep='-')
  # add.to.list (files, annotation)
}

##
## run similarity calcs and visualize matrix
##
similarity (master.list, 'overall-similarity-matrix', master.annotation)


