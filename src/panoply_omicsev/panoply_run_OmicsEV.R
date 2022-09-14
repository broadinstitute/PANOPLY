library(OmicsEV)

################################################################################
# handle command line arguments
cat("\nextracting command line arguments\n")

args = commandArgs(trailingOnly = T)

if (length(args) == 6) {
  
  # all inputs provided
  data_dir <- args[1]
  sample_list <- args[2]
  cpu <- as.integer(args[3])
  data_type <- args[4]
  x2 <- args[5]
  do_fun_pred <- args[6]
  
} else {
  stop("Incorrect number of arguments")
}

if (!file.exists(x2)) {
  x2 <- NULL # no rna file for correlation
}

################################################################################

cat("data_dir:", data_dir, '\n')
cat("sample_list:", sample_list, '\n')
cat("x2:", x2, '\n')
cat("cpu:", cpu, '\n')
cat("data_type:", data_type, '\n')
cat('do_fun_pred:', do_fun_pred, '\n\n')

################################################################################

print("running the full function")

run_omics_evaluation(data_dir = data_dir,
                     sample_list = sample_list,
                     x2 = x2,
                     cpu=cpu,
                     data_type= data_type,
                     do_fun_pred = do_fun_pred)
