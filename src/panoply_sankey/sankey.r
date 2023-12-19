#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
### Create Rmarkdown report for consensus clustering module ###

# tar_file  - URL of tar file created by task panoply_cons_clust
# yaml_file - URL of master parameters yaml file including output from startup notebook
# label     - character, name of folder in tarball
# data_type - character, data type

args = commandArgs(TRUE)

label = args[1]

annot_files_str = args[2] # filenames, concatenated into a string, separated by commas (e.g. "file_a.txt,file_b.txt")
annot_file_types_str = args[3] # filetype labels, concatenated into a string, separated by commas (e.g. "Proteome,Phosphoproteome")

annot_column = args[4] # annotation column of interest (e.g. "NMF.consensus")
annot_prefix = args[5] # prefix that should be appended to annotation-entries (e.g. 'C' for clusters-- C1, C2, C3 instead of 1, 2, 3)
if (is.na(annot_prefix)) annot_prefix=""

sample_id_col = "Sample.ID"

annot_file_primary = args[6] # optional "primary comparison" file (e.g. Multiomic data), which should be highlighted / centered in comparisons
annot_file_type_primary = args[7] # optional; filetype label for "primary" comparison


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

# For NMF-related stuff (I think)
source('https://raw.githubusercontent.com/karstenkrug/R-code/main/my_plots.r')
#source('c:/Users/karsten/Dropbox/Devel/PANOPLY/src/panoply_mo_nmf/nmf_functions.R')
source('https://raw.githubusercontent.com/broadinstitute/PANOPLY/dev/src/panoply_mo_nmf/nmf_functions.R')




### Begin Editable Section -----

annot_files = unlist(strsplit(annot_files_str, ","))
annot_file_types = unlist(strsplit(annot_file_types_str, ","))

# sanity check that we have the same number of files and labels
if (length(annot_files)!=length(annot_file_types)) stop( glue("The number of files ({length(annot_files)}) provided for comparison does not match the number of labels ({length(annot_file_types)}) provided."))

names(annot_files) = annot_file_types


if (file.exists(annot_file_primary)) { #if we have a secondary datafile to compare to
  annot_file_types = c(annot_file_type_primary, annot_file_types) #add label to array (at BEGINNING)
  annot_files[annot_file_type_primary] = annot_file_primary #add file to array
}


# combos to plot (all combos either order)
datatype_combos <- c(combn(annot_file_types, 2, simplify = FALSE), #forwards
                combn(annot_file_types, 2, rev, simplify = FALSE)) #backwards

if (file.exists(annot_file_primary) && length(annot_file_types)>2) { #if we have a secondary datafile to compare to, and have at least three files total
  datatype_combos = c(datatype_combos,
                 combn(annot_file_types[which(annot_file_types!=annot_file_type_primary)], 2, #take all datatype combos besides annot_file_type_primary
                       FUN = function(arr) {append(arr, annot_file_type_primary, after=1)} , simplify = FALSE)) #add annot_file_type_primary in between eaach combo
}


### End Editable Section -----




###################################
## import files and select columns
data_list <- lapply(annot_files, function(filename) { # read in / validate
  data = read_delim(filename, delim='\t')
  if (!(sample_id_col %in% colnames(data))) stop(paste(glue("Could not locate sample_id_col ('{sample_id_col}')."))) # check for sample ID column
  if (!(annot_column %in% colnames(data))) stop(paste(glue("Could not locate annot_column ('{annot_column}')."))) # check for annotation column
  return(data)
})
data_subset_list <- lapply(data_list, function(data) { # subset data
  data_subset = data %>% 
    select( all_of( c( sample_id_col, annot_column))) %>% # select sample id
    mutate(!!sample_id_col := make.unique(!!sym(sample_id_col))) %>% # force unique sample IDs
    mutate(!!annot_column := paste0(annot_prefix, !!sym(annot_column)) ) %>% # add prefix to annot_column, if desired
    mutate(!!annot_column := factor(!!sym(annot_column), levels=sort(unique(!!sym(annot_column))))) # factorize annot_column, 
  return(data_subset)
}) 


##################################
## set file_type as column name for annot_column
data_subset_list <- lapply(names(data_subset_list), function(data_label){
  rename(data_subset_list[[data_label]],
         !!data_label := annot_column) # rename column with data_label
} )

#######################################
## merge tables
my_join <- function(x, y){
  full_join(x, y, by=sample_id_col)
  #inner_join(x, y, by=sample_id_col)
}
data_full <- Reduce(my_join, data_subset_list) %>% # merge dataframes by sample ID
  column_to_rownames(., sample_id_col) %>% data.frame # set sample IDs to rownames

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
tar(paste(label,'sankey_diagrams.tar.gz', sep="_"), list.files(pattern="sankey-.+?.html"))
