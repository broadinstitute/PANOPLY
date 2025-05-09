library(mice)
library(cmapR)
library(optparse)
library(dplyr)
library(glue)

## read command line arguments
option_list <- list(
  make_option(c('-f', '--file_path'), type = 'character'),
  make_option(c('-n', '--na_max'), type = 'numeric'),
  make_option(c('-m', '--num_imputations'), type = 'integer', default = 15L),
  make_option(c('-c', '--num_cores'), type = 'integer', dest = 'num_cores', default = 1L),
  make_option(c('-s', '--seed'), type = 'integer', dest = 'seed', default = 1L),
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.')
  # make_option(c('-r', '--rdata_out_file_name'), type = 'character', default = 'imputation_mice_object.RData')
)

opt <- parse_args(OptionParser(option_list = option_list))
file_path <- opt$file_path
na_max <- opt$na_max
num_imputations <- opt$num_imputations
num_cores <- opt$num_cores
seed <- opt$seed
output_prefix <- opt$output_prefix
# out_file_name <- opt$out_file_name
# rdata_out_file_name <- opt$rdata_out_file_name




#### Natalie's Version ####
#impute missing values using mice
mice_imputation <- function(gct_file, na_max=0.4, num_imps=15, seed=2023, num_cores=1, output_prefix = "results"){
  
  print(paste("## na_max:", na_max))
  print(paste("## num_imps:", num_imps))
  print(paste("## seed:", seed))
  print(paste("## num_cores:", num_cores))
  
  #code from Stephanie Vartany to run Mice
  # read in data (features x samples)
  # typically MICE does samples x features (feature-wise), but this takes much too long. Stephanie has had good results with sample-wise which we perform here. With parallelization, feature-wise may be possible, but is likely not worth it.
  gct <- parse_gctx(gct_file)
  data <- gct@mat
  
  #filter using na_max
  #MICE seems to have optimal results when na_max=0.4, after that it starts to drop off
  #this is likely dataset dependent
  keep <- row.names(data)[rowSums(is.na(data))/dim(data)[2] <= na_max]
  gct_filt <- subset_gct(gct,rid=keep)
  data_filt <- gct_filt@mat
  
  print(glue("## Filtered GCT to features with <{round(na_max*100,2)}% missing values. {dim(data)[1]-dim(data_filt)[1]} features dropped for incompleteness."))
  
  
  # run mice, this creates a mice-specific output object
  # m is number of iterations, set to 15 (default was 5)
  # in mice v3.15 and later you can run in parallel using futuremice()
  # set seed so result is the same for the same dataset (for futuremice, use parallelseed)
  print("Imputing using MICE")
  if(num_cores>1){
    mice_out <- futuremice(data_filt, m=num_imps, parallelseed=seed, n.core=num_cores, print=T)
  }else{
    mice_out <- mice(data_filt, m=num_imps, seed=seed, print=T)
  }
  print("Imputation done")
  
  # collect the mice object (with all iterations)
  print("Aggregating imputations")
  imp_data <- mice::complete(mice_out, 'long')
  
  # aggregate across all iterations
  avg_imp_data <- imp_data %>% 
    group_by(.id) %>% 
    summarize(across(.fns = mean)) %>%
    select(-c('.id', '.imp'))
  avg_imp_data <- as.data.frame(avg_imp_data)
  rownames(avg_imp_data) <- rownames(data_filt)
  print("Aggregation done")
  
  #save avg_imp_data to .csv file
  print("Saving imputed data")
  fn <- paste(output_prefix,"_mice_imputed.csv",sep="")
  write.csv(avg_imp_data,fn)
  
  #save workspace
  fn <- paste(output_prefix,"_mice_imputed.RData",sep="")
  save.image(fn)
  
  # save avg_imp_data to gct file
  gct@mat <- as.matrix(avg_imp_data)
  fn <- paste(output_prefix,"_mice_imputed.gct",sep="")
  write_gct(gct,fn,appenddim = F)
  print("Saving done")
  
}

mice_imputation(gct_file = file_path,
                na_max = na_max,
                num_imps = num_imputations,
                seed = seed,
                num_cores = num_cores,
                output_prefix = output_prefix)


# #### Stephanie's Version ####
# 
# 
# # read in file
# ext <- tools::file_ext(file_path)
# if (ext == 'csv') {
#   data <- read.csv(file_path, row.names = 1)
# } else if (ext == 'gct') {
#   data <- as.data.frame(parse_gctx(file_path)@mat)
# }
# 
# # check if data has any missing values
# if (!any(is.na(data))) {
#   warning("Data does not contain any missing values, MICE will not be run.")
#   avg_imp_data <- data
#   mice_out <- NULL
#   
# } else {
#   # run mice, this creates a mice-specific output object
#   # m is number of iterations
#   cat(paste("\nRunning MICE with", num_imputations, "imputations.\n"))
#   mice_out <- futuremice(data,
#                          m = num_imputations,
#                          n.core = num_cores,
#                          parallelseed = seed)
#   
#   # collect the mice object (with all iterations)
#   imp_data <- mice::complete(mice_out, 'long')
#   
#   # aggregate across all iterations
#   avg_imp_data <- imp_data %>%
#     filter(.imp <= num_imputations) %>%
#     group_by(.id) %>%
#     summarize(across(.cols = everything(), .fns = mean)) %>%
#     select(-c('.id', '.imp'))
#   avg_imp_data <- as.data.frame(avg_imp_data)
#   rownames(avg_imp_data) <- rownames(data)
#   
#   # avg_imp_data should be your final output
#   cat("\nMICE imputation complete!\n")
# }
# 
# cat("\nSaving outputs...")
# 
# # save the outputs
# write.csv(avg_imp_data, file = out_file_name, row.names = T)
# save(list = 'mice_out', file = rdata_out_file_name)
# 
# cat("DONE\n\n")
