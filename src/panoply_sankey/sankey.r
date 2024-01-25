#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
### Create Sankey Diagrams ###

args = commandArgs(TRUE)


rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  #### NMF Outputs ####
  make_option( c("-f", "--annot_files"), action='store', type='character', dest='annot_files_str', help='String of TSV annotation-files, separated by commas.') , 
  make_option( c("-l", "--annot_file_labels"), action='store', type='character', dest='annot_file_types_str', help='String of labels for each annotation-file, separated by commas.') , 
  make_option( c("-j", "--annot_file_primary"), action='store', type='character', dest='annot_file_primary', help='TSV annotation-file which should be centered / highlighted in comparisons.') ,
  make_option( c("-m", "--annot_label_primary"), action='store', type='character', dest='annot_file_type_primary', help='Labels for primary annotation file.') , 
  #### Annotation File Format ####
  make_option( c("-i", "--id_column"), action='store', type='character',  dest='id_column', help='ID column which uniquely identifies entries. Must be shared across all annotation files.'),
  make_option( c("-a", "--annot_column"), action='store', type='character',  dest='annot_column', help='Column annotation to create sankey comparisons for.'),
  make_option( c("-p", "--annot_prefix"), action='store', type='character',  dest='annot_prefix', default="", help='Prefix to prepend to annotations (e.g. \'C\' for C1, C2, C3 instead of 1, 2, 3).'),
  #### General Parameters ####
  make_option( c("-x", "--label"), action='store', type='character',  dest='label', help='Label associated with this run.')
  #### ####
)

opt <- parse_args( OptionParser(option_list=option_list),
                   # # for testing arguments
                   # args = c('--annot_files',"/opt/input/odg_all-so_nmf-prot_K3_clusterMembership.tsv,/opt/input/odg_all-so_nmf-pSTY_K3_clusterMembership.tsv,/opt/input/odg_all-so_nmf-acK_K4_clusterMembership.tsv,/opt/input/odg_all-so_nmf-ubK_K5_clusterMembership.tsv",
                   #          '--annot_file_labels',"prot,pSTY,acK,ubK",
                   #          '--annot_file_primary',"/opt/input/odg_all-mo_nmf_K4_clusterMembership.tsv",
                   #          '--annot_label_primary',"Multiomic",
                   #          '--annot_column',"NMF.consensus",
                   #          '-x',"odg_test")
)


library(pacman)
p_load(readr)
p_load(dplyr)
p_load(purrr)
p_load(glue)
p_load(tibble)

library(htmlwidgets) #must be 1.2.0, don't try to p_load
# p_load(webshot)
# webshot::instfull_phantomjs()


# For SankeyDiagram()
library('flipU')
library('flipTransformations')
p_load('networkD3') # needed for sankeyNetwork()
source('https://raw.githubusercontent.com/Displayr/flipPlots/master/R/sankeydiagram.R')
# NOTE: the full flipPlots package requires a newer version of R, which is why ONLY the sankeydiagram.R file is loaded

# # For NMF-related stuff (I think)
# source('https://raw.githubusercontent.com/karstenkrug/R-code/main/my_plots.r')
# #source('c:/Users/karsten/Dropbox/Devel/PANOPLY/src/panoply_mo_nmf/nmf_functions.R')
# source('https://raw.githubusercontent.com/broadinstitute/PANOPLY/dev/src/panoply_mo_nmf/nmf_functions.R')




### Begin Editable Section -----

annot_files = unlist(strsplit(opt$annot_files_str, ","))
annot_file_types = unlist(strsplit(opt$annot_file_types_str, ","))

# sanity check that we have the same number of files and labels
if (length(annot_files)!=length(annot_file_types)) stop( glue("The number of files ({length(annot_files)}) provided for comparison does not match the number of labels ({length(annot_file_types)}) provided."))

names(annot_files) = annot_file_types


if (!is.null(opt$annot_file_primary)) { #if we have a secondary datafile to compare to
  annot_file_types = c(opt$annot_file_type_primary, annot_file_types) #add label to array (at BEGINNING)
  annot_files[opt$annot_file_type_primary] = opt$annot_file_primary #add file to array
}


# combos to plot (all combos either order)
datatype_combos <- c(combn(annot_file_types, 2, simplify = FALSE), #forwards
                     combn(annot_file_types, 2, rev, simplify = FALSE)) #backwards

if (!is.null(opt$annot_file_primary) && length(annot_file_types)>2) { #if we have a secondary datafile to compare to, and have at least three files total
  datatype_combos = c(datatype_combos,
                 combn(annot_file_types[which(annot_file_types!=opt$annot_file_type_primary)], 2, #take all datatype combos besides annot_file_type_primary
                       FUN = function(arr) {append(arr, opt$annot_file_type_primary, after=1)} , simplify = FALSE)) #add annot_file_type_primary in between eaach combo
}


### End Editable Section -----




###################################
## import files and select columns
data_list <- lapply(annot_files, function(filename) { # read in / validate
  cat(glue("Reading in '{filename}'.\n\n"))
  data = read.delim(filename)
  if (is.null(opt$id_column)) {  # if we have no ID column
    if (! all(rownames(data)==1:length(rownames(data))) ) { # unless we have default numeric rownames
      opt$id_column = "sankey_id_column" # create an ID column
      data[[opt$id_column]] = rownames(data) # use rownames as ID column
    } else {
      stop(paste(glue("No ID column selected, and dataset does not contain appropriate rownames')."))) # check for sample ID column
    }
  } # if we have no ID column, use rownames or complain
  if (!(opt$id_column %in% colnames(data))) stop(paste(glue("Could not locate sample_id_col ('{opt$id_column}')."))) # check for sample ID column
  if ( sum(duplicated(data[[opt$id_column]]))>=1 ) stop(paste(glue("The identified column contains duplicated IDs. Please provide an id_column which uniquely identifies entries."))) # check for sample ID column
  if (!(opt$annot_column %in% colnames(data))) stop(paste(glue("Could not locate annot_column ('{opt$annot_column}')."))) # check for annotation column
  return(data)
})
# set ID column, based on whether we're using rownames, or it was set with opt$id_column  
if (is.null(opt$id_column)) {
  id_column = "sankey_id_column" # if the ID column wasn't specified, overwrite with the generated ID column
} else {
  id_column = opt$id_column
}
# subset data to columns of interest / add prefix if specified
data_subset_list <- lapply(data_list, function(data) { # subset data
  data_subset = data %>% 
    select( all_of( c( id_column, opt$annot_column))) %>% # select sample id
    # mutate(!!id_column := make.unique(!!sym(id_column))) %>% # force unique sample IDs
    mutate(!!opt$annot_column := paste0(opt$annot_prefix, !!sym(opt$annot_column)) ) %>% # add prefix to annot_column, if desired
    mutate(!!opt$annot_column := factor(!!sym(opt$annot_column), levels=sort(unique(!!sym(opt$annot_column))))) # factorize annot_column, 
  return(data_subset)
}) 


##################################
## set file_type as column name for annot_column
data_subset_list <- lapply(names(data_subset_list), function(data_label){
  rename(data_subset_list[[data_label]],
         !!data_label := opt$annot_column) # rename column with data_label
} )

#######################################
## merge tables
my_join <- function(x, y){
  full_join(x, y, by=id_column)
  #inner_join(x, y, by=id_column)
}
data_full <- Reduce(my_join, data_subset_list) %>% # merge dataframes by sample ID
  column_to_rownames(., id_column) %>% data.frame # set sample IDs to rownames

# this appends the column name to the cluster (e.g.  acK_C1 vs C1) and is very redundant
# for(i in 1:ncol(data_full))
#     data_full[, i] <- paste(colnames(data_full)[i], '_', data_full[, i], sep='')


#######################################
# reformat data into sankey-interpretable input
cont.tab.flat <- apply(data_full, 1, paste, collapse='-') %>% # 'flatten' datatable into connections
  table(.) # sum connections

cont.tab.split <- strsplit(names(cont.tab.flat), '-') # pull out ONLY the connections that actually exist
#cont.tab.clust <- lapply(cont.tab.split, function(x) sub('.*_','', x))
cont.tab.clust <- cont.tab.split

my.data <- data.frame(Reduce('rbind', cont.tab.clust), as.numeric(cont.tab.flat)) # convert to data.frame
colnames(my.data) <- c(names(data_full), 'freq') # add filetype labels back into data

max_annots=sum(unique(unlist(my.data[-length(my.data)]))!='NA') #figure out the max number of annots, across all files, to enforce consistent colors

for (datatypes_of_interest in datatype_combos) {
  
  colors=colorRampPalette(c('#fde0dd','#f768a1','#7a0177'))(max_annots) # choose enough colors to color all clusters

  ###################################
  ## plot
  # link.color = 'Source'
  for (link.color in c('Source', 'Target')) {
    widget <- SankeyDiagram(my.data[, datatypes_of_interest],
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
    datatype_tmp <- paste(datatypes_of_interest, collapse='_')
    fn <- glue("sankey-{datatype_tmp}-N-{sum(my.data$freq)}-{link.color}.html")
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
tar(paste(opt$label,'sankey_diagrams.tar.gz', sep="_"), list.files(pattern="sankey-.+?.html"))
