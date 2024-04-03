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
  make_option( c("-c", "--color_str"), action='store', type='character',  dest='color_str', help='String of hex-codes for sankey diagram colors, separated by commas. If too many/few colors are provided, colors will be extrapolated.'),
  make_option( c("-x", "--label"), action='store', type='character',  dest='label', help='Label associated with this run.')
  #### ####
)

opt <- parse_args( OptionParser(option_list=option_list),
                   # # for testing arguments
                   # args = c('--annot_files',"/opt/input/melanoma-v1-so_nmf-prot_K3_clusterMembership.tsv,/opt/input/melanoma_v1_DIA-mo_nmf_K4_clusterMembership.tsv,/opt/input/melanoma_v1_k=4_250iter-mo_nmf_K4_clusterMembership.tsv,/opt/input/melanoma_v1_DIA-so_nmf-DIA_prot_K4_clusterMembership.tsv,/opt/input/melanoma-v1-so_nmf-prot_K3_clusterMembership.tsv,/opt/input/melanoma_v1_DIA-so_nmf-RNA_K3_clusterMembership.tsv",
                   #          '--annot_file_labels',"Multi_DIA_k=3,Multi_DIA_k=4,Multi_DDA,DIA_prot,DDA_prot,RNA" ,
                   #          '--annot_file_primary',"/opt/input/melanoma-v1-so_nmf-prot_K3_clusterMembership.tsv",
                   #          '--annot_label_primary',"Multi_DIA_k=3",
                   #          '--annot_column',"NMF.consensus",
                   #          '-c',"#AA0000,#0000AA",
                   #          '-x',"melanoma_v1_DIA_vs_DDA")
)


# library(pacman)
library(readr)
library(dplyr)
library(purrr)
library(glue)
library(tibble)



# For SankeyDiagram()
library('flipU')
library('flipTransformations')
library('networkD3') # needed for sankeyNetwork()
source('https://raw.githubusercontent.com/Displayr/flipPlots/master/R/sankeydiagram.R')
# NOTE: the full flipPlots package requires a newer version of R, which is why ONLY the sankeydiagram.R file is loaded

library(htmlwidgets) #must be 1.2.0, don't try to p_load
library(webshot)


# # For NMF-related stuff (I think)
# source('https://raw.githubusercontent.com/karstenkrug/R-code/main/my_plots.r')
# #source('c:/Users/karsten/Dropbox/Devel/PANOPLY/src/panoply_mo_nmf/nmf_functions.R')
# source('https://raw.githubusercontent.com/broadinstitute/PANOPLY/dev/src/panoply_mo_nmf/nmf_functions.R')




###################################
## specify annotation files and their respective labels

annot_files = unlist(strsplit(opt$annot_files_str, ",")) # split into array of filenames
annot_file_types = unlist(strsplit(opt$annot_file_types_str, ",")) # split into array of labels #%>% make.names(unique = TRUE) # remove special characters

# sanity check that we have the same number of files and labels
if (length(annot_files)!=length(annot_file_types)) stop( glue("The number of files ({length(annot_files)}) provided for comparison does not match the number of labels ({length(annot_file_types)}) provided."))

names(annot_files) = annot_file_types

# establish a primary annotation file, if provided
if (!is.null(opt$annot_file_primary)) { #if we have a secondary datafile to compare to
  annot_file_types = c(opt$annot_file_type_primary, annot_file_types) #%>% #add label to array (at BEGINNING)
    #make.unique() # ensure that labels are unique
  annot_files[opt$annot_file_type_primary] = opt$annot_file_primary #add file to array
}
# save labels to file, to avoid regex parsing in the report
readr::write_lines(annot_file_types, 'sankey_labels.txt')



###################################
## specify annotation-combinations we want to create sankey diagrams for
datatype_combos <- c(combn(annot_file_types, 2, simplify = FALSE), #forwards
                     combn(annot_file_types, 2, rev, simplify = FALSE)) #backwards

if (!is.null(opt$annot_file_primary) && length(annot_file_types)>2) { #if we have a secondary datafile to compare to, and have at least three files total
  datatype_combos = c(datatype_combos,
                 combn(annot_file_types[-1], 2, #take all datatype combos besides annot_file_type_primary (i.e. "first" filetype)
                       FUN = function(arr) {append(arr, opt$annot_file_type_primary, after=1)} , simplify = FALSE)) #add annot_file_type_primary in between eaach combo
}




###################################
## import files and select columns
data_list <- lapply(annot_files, function(filename) { # read in / validate
  cat(glue("Reading in '{filename}'.\n\n"))
  if (is.null(opt$id_column)) {  # if we have no ID column
    data = read.delim(filename, row.names = 1) # assume rownames
    if (! all(rownames(data)==1:length(rownames(data))) ) { # unless we have default numeric rownames
      opt$id_column = "sankey_id_column" # create an ID column
      data[[opt$id_column]] = rownames(data) # use rownames as ID column
    } else {
      stop(paste(glue("No ID column selected, and dataset does not contain appropriate rownames')."))) # check for sample ID column
    }
  } else {
    data = read.delim(filename)
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
  column_to_rownames(., id_column)# %>% data.frame # set sample IDs to rownames


#######################################
# reformat data into sankey-interpretable input
cont.tab.flat <- apply(data_full, 1, paste, collapse='-') %>% # 'flatten' datatable into connections
  table(.) # sum connections

cont.tab.split <- strsplit(names(cont.tab.flat), '-') # pull out ONLY the connections that actually exist
#cont.tab.clust <- lapply(cont.tab.split, function(x) sub('.*_','', x))
cont.tab.clust <- cont.tab.split

my.data <- data.frame(Reduce('rbind', cont.tab.clust), as.numeric(cont.tab.flat)) # convert to data.frame
colnames(my.data) <- c(names(data_full), 'freq') # add filetype labels back into data

#######################################
# set color-scheme for annotations
# get all annotation values
annots = unique(unlist(my.data[-length(my.data)]))
annots_noNA = annots[annots!='NA'] # exclude NA
max_annots=length(annots_noNA) # count number of non-NA values
# format color-scheme
default_palette = c('#fde0dd','#f768a1','#7a0177') # default color-palette
if(!is.null(opt$color_str)) { # if the user provided colors as an argument
  palette = unlist(strsplit(opt$color_str, ",")) # parse the argument string into a vector
} else {
  palette = default_palette # otherwise use the default palette
}
# create color-scale for annot values
colors = tryCatch(colorRampPalette(palette)(max_annots),
                  error = function(e) {
                    cat(glue("\nWARNING: Some or all of the provided colors ({paste(palette, collapse=', ')}) were invalid. Using default palette.\n"))
                    colorRampPalette(default_palette)(max_annots) # choose enough colors to color all clusters
                  })
# names(colors) = annots_noNA # name colors (SankeyDiagram doesn't consider names)
# colors = c(colors, "NA" = "grey") # add NA color explicitly (SankeyDiagram doesn't consider names)

for (datatypes_of_interest in datatype_combos) {
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
    # save widget to HTML file
    saveWidget(widget, file=fn, selfcontained = TRUE, libdir = NULL,
               background = "white", knitrOptions = list())
    # # save widget to PDF file
    webshot(fn, file=sub('\\.html$','.pdf',fn)) # takes a screenshot of the sankey diagram -> PDF
    # cmd <- sprintf("pandoc '%s' --pdf-engine=pdflatex -t latex -o '%s'", fn, sub('\\.html$','.pdf',fn)) # create pandoc command-line prompt for generating PDF
    # system(cmd) # run command-line prompt to convert HTML to PDF 
    
    file.copy(dir('.', pattern = fn), wd, overwrite = T)
    file.copy(dir('.', pattern = sub('\\.html$','.pdf',fn)), wd, overwrite = T)
    setwd(wd)
  }
  
}

# # tar sankey_diagrams into a single file
# out_tar = paste(opt$label,'sankey_diagrams.tar.gz', sep="_") # tar file-name
# out_files = c(list.files(pattern="sankey-.+?(.html)|(.pdf)"), # get all HTML/PDFs filenames
#               "sankey_labels.txt") # and also get the TXT with the file-labels, so the report doesn't need regex
# cmd = glue("tar -czvf {out_tar} {paste(out_files, collapse=' ')}") # run tar as a command-line prompt, since tar() is failing
# system(cmd) # run command-line tar prompt


# tar(tarfile = ,
#     files = list.files(pattern="sankey-.+?(.html)|(.pdf)"), # collect all HTML and PDF files
#     compression = "gzip", # compress using gzip (i.e. .gz file)
#     extra_flags = "-v") # verbose
