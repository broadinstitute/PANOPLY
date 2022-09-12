library(OmicsEV)

################################################################################
# handle command line arguments
cat("extracting command line arguments")

args = commandArgs(trailingOnly = T)

if (length(args) == 5) {
  
  # all inputs provided
  data_dir <- args[1]
  sample_list <- args[2]
  cpu <- as.integer(args[3])
  data_type <- args[4]
  x2 <- args[5]
  
} else {
  stop("Incorrect number of arguments")
}

if (!file.exists(x2)) {
  x2 <- NULL # no rna file for correlation
}

################################################################################

cat(paste("data_dir:", data_dir, sep=' '))
cat(paste("sample_list:", sample_list, sep=' '))
cat(paste("x2:", x2, sep=' '))
cat(paste("cpu:", cpu, sep=' '))
cat(paste("data_type:", data_type, sep=' '))

################################################################################

print("running the full function")

run_omics_evaluation(data_dir = data_dir,
                     sample_list = sample_list,
                     x2 = x2,
                     cpu=cpu,
                     data_type= data_type)
