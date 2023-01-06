#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
### Create Rmarkdown report for consensus clustering module ###

# tar_file  - URL of tar file created by task panoply_cons_clust
# yaml_file - URL of master parameters yaml file including output from startup notebook
# label     - character, name of folder in tarball
# data_type - character, data type

args = commandArgs(TRUE)

so_tar_file = args[1]
label = args[2]
mo_tar_file = args[3]

library(pacman)
p_load(readr)
p_load(dplyr)
p_load(purrr)
p_load(glue)
p_load(tibble)

p_load(htmlwidgets)
# p_load(webshot)
# webshot::instfull_phantomjs()


# For SankeyDiagram()
library('flipTransformations')
p_load('networkD3') # needed for sankeyNetwork()
source('https://raw.githubusercontent.com/Displayr/flipPlots/master/R/sankeydiagram.R')
# NOTE: the full flipPlots package requires a newer version of R, which is why ONLY the sankeydiagram.R file is loaded

# For NMF-related stuff (I think)
source('https://raw.githubusercontent.com/karstenkrug/R-code/main/my_plots.r')
#source('c:/Users/karsten/Dropbox/Devel/PANOPLY/src/panoply_mo_nmf/nmf_functions.R')
source('https://raw.githubusercontent.com/broadinstitute/PANOPLY/dev/src/panoply_mo_nmf/nmf_functions.R')




### Begin Editable Section -----
# please note that the sample ID column is required to be Sample.ID

# provide path for each data set's clin_annot_nmf.txt file
# variable naming scheme isn't important
# these can be found in the tar outputted by NMF clustering

### Unpack SO-NMF results

untar(so_tar_file)
data.dir = 'nmf_results/so-nmf' # directory structure within the tar, from panoply_so_nmf_assemble 

ome_dirs = list.files(data.dir, full.names = TRUE)
ome_types = stringr::str_extract(ome_dirs, '([^-]+?)$') #extract ome-type from directory names
names(ome_dirs) = ome_types #name ome_dirs with ome_types

### Location of clin_annot_nmf.txt files for SO analysis
ome_clust = c()
for (ome in ome_types) {
  ome_clust[ome] = list.files(ome_dirs[ome], pattern = 'clin_anno_nmf.txt', full.names = TRUE, recursive=TRUE)
}


if (file.exists(mo_tar_file)) { #if we have Multiomic data
  ome_types = c("Multiomic", ome_types) #add Multiomic as an ome-type
  
  untar(mo_tar_file, exdir="mo_nmf")
  mo_file = list.files("mo_nmf", pattern='clin_anno_nmf.txt', full.names =TRUE, recursive=TRUE) #locate clin_annot file

  ome_clust["Multiomic"] = mo_file #add file to ome_clust files
}


# combos to plot (all combos either order)
ome_combos <- c(combn(ome_types, 2, simplify = FALSE), #forwards
                combn(ome_types, 2, rev, simplify = FALSE)) #backwards

if (file.exists(mo_tar_file)) { #if we have Multiomic data
  ome_combos = c(ome_combos,
                 combn(ome_types[which(ome_types!="Multiomic")], 2, #take all ome combos besides Multiomic
                       FUN = function(arr) {append(arr, "Multiomic", after=1)} , simplify = FALSE)) #add Multiomic in between eaach combo
}


### End Editable Section -----




###################################
## import
tab <-  lapply(ome_clust, read_delim, delim='\t') %>%
  lapply(., function(x){
    x$NMF.consensus = paste0('C', x$NMF.consensus)
    x$NMF.consensus.core[!is.na(x$NMF.consensus.core)] <- paste0('C', x$NMF.consensus.core[!is.na(x$NMF.consensus.core)])
    #x$NMF.consensus.core[is.na(x$NMF.consensus.core)] <- paste0(x$NMF.consensus[is.na(x$NMF.consensus.core)], '-mixed')
    return(x)
  })

tab <- lapply(tab, function(x) {
  if (!('Sample.ID' %in% colnames(x))) #if we don't have Sample.ID, 
    { warning(paste(glue("Could not locate Sample.ID column; assuming the first column has Sample.IDs.")))
      colnames(x)[1]='Sample.ID'; return(x)} else # print a WARNING and assume the first column has Sample.ID
    { colnames(x)[1]='idddd'; return(x)} }) #otherwise, name it something that won't overlap anything
tab <- lapply(tab, function(x)
  x %>% 
    #filter(!grepl('mixed', NMF.consensus.core)) %>%
    select( one_of( c( 'Sample.ID', 'NMF.consensus'))) %>%
    mutate(Sample.ID=make.unique(Sample.ID))
) 

#tab_tmp <- tab
##for(i in names(tab_tmp)){
#  tab_tmp[[i]]$NMF.consensus <- paste( tab_tmp[[i]]$NMF.consensus, i, sep='_')
#}

#tab[[1]] <- Reduce('rbind', tab_tmp)
#names(tab)[1] <- 'conf-disc'
#tab <- tab[-c(2)]


##################################
## append data type to column names
tab <- lapply(names(tab), function(x){
  xx <- tab[[x]]
  colnames(xx)[2] <- glue("{x}")
  xx
} )

#######################################
## join
my_join <- function(x, y){
  full_join(x, y, by='Sample.ID')
  #inner_join(x, y, by='Sample.ID')
}

tab <- Reduce(my_join, tab)
tab <- column_to_rownames(tab, 'Sample.ID') %>% data.frame

# this appends the column name to the cluster (e.g.  acK_C1 vs C1) and is very redundant
# for(i in 1:ncol(tab))
#     tab[, i] <- paste(colnames(tab)[i], '_', tab[, i], sep='')

cont.tab.flat <- apply(tab, 1, paste, collapse='-')
cont.tab.flat <- table(cont.tab.flat)


cont.tab.split <- strsplit(names(cont.tab.flat), '-')
#cont.tab.clust <- lapply(cont.tab.split, function(x) sub('.*_','', x))
cont.tab.clust <- cont.tab.split

my.data <- data.frame(Reduce('rbind', cont.tab.clust), as.numeric(cont.tab.flat))
colnames(my.data) <- c(names(tab), 'freq')

max_nclust=sum(unique(unlist(my.data[-length(my.data)]))!='NA') #figure out the max number of clusters
my.data[-length(my.data)] = lapply(my.data[-length(my.data)], factor, levels=paste0('C',1:max_nclust)) #factor clusters so they're ordered properly

for (omes_of_interest in ome_combos) {
  
  colors=colorRampPalette(c('#fde0dd','#f768a1','#7a0177'))(max_nclust) # choose enough colors to color all clusters

  ###################################
  ## plot
  # link.color = 'Source'
  for (link.color in c('Source', 'Target')) {
    widget <- SankeyDiagram(my.data[, omes_of_interest],
                            link.color = link.color,
                            # link.color = "First variable",
                            # link.color = "None",
                            
                            weights = my.data$freq,
                            variables.share.values = TRUE,
                            colors=colors,
                            # label.show.percentages = TRUE,
                            node.padding=50,font.size = 20,
                            max.categories = 15)
    
    # widget    
    
    # ## export
    ome_tmp <- paste(omes_of_interest, collapse='_')
    fn <- glue("sankey-{ome_tmp}-N-{sum(my.data$freq)}-{link.color}.html")
    tmp.dir <- tempdir()
    
    wd <- getwd()
    setwd(tmp.dir)
    saveWidget(widget, file=fn, selfcontained = TRUE, libdir = NULL,
               background = "white", knitrOptions = list())
    # webshot(fn, file=sub('\\.html$','.pdf',fn)) #would take a screenshot of the sankey diagram -> PDF. don't need
    file.copy(dir('.', pattern = fn), wd, overwrite = T)
    file.copy(dir('.', pattern = sub('\\.html$','.pdf',fn)), wd, overwrite = T)
    setwd(wd)
  }
  
}

# tar sankey_diagrams into a single file
tar(paste(label,'sankey_diagrams.tar.gz', sep="_"), list.files(pattern="sankey-.+?.html"))