library(dplyr)

combine_methods=function(method1_folder, method2_folder, 
                         sample_annotation_file,
                         clinical_attributes = c("gender"), 
                         out_dir = "./", prefix = "cosmo"){
  
  #sj_intermediate_file <- paste(method1_folder,"/intermediate.csv", sep = "")
  sj_intermediate_file <- paste(method1_folder,"/clinical_attributes_pred.tsv", sep = "")
  sj_correlation_file <- paste(method1_folder,"/sample_correlation.csv", sep = "")
  modelA_result_file <- paste(method2_folder,"/test_ModelA_result.csv", sep = "")
  modelB_result_file <- paste(method2_folder,"/test_ModelB_result.csv", sep = "")
  
  cli_data_file <- sample_annotation_file
  cli_data <- read.delim(cli_data_file,stringsAsFactors = FALSE)
  
  sention_d1 <- read.table(modelA_result_file, sep=',', header=TRUE)
  sention_d2 <- read.table(modelB_result_file, sep=',', header=TRUE)
  
  sjcli <- read.table(sj_intermediate_file, header=TRUE, sep = "\t")
  corsample <- read.table(sj_correlation_file, sep=',', header=TRUE, row.names=1)
  corsample <- as.matrix(corsample)
  
  
  rankdist <- computeRankdist(corsample)
  malerank <- rankbyRow(rankdist)
  fmlerank <- rankbyCol(rankdist)
  
  matcher <- stableMarriage(malerank, fmlerank)
  d1match <- c()
  d2match <- c()
  for(i in 1:length(matcher$malematch)){
    d1match[i] <- which(rownames(corsample) == matcher$malematch[i])
  }
  
  for(i in 1:length(matcher$fmlematch)){
    d2match[i] <- which(rownames(corsample) == matcher$fmlematch[i])
  }
  
  
  total_samples <- nrow(sjcli)
  nonmatch <- which(1:total_samples != d1match)
  length(nonmatch)
  
  
  ### get pairdist for shifting chain identification later
  pairdist <- data.frame(d1=1:total_samples, d2=d1match)
  pairdist$d1rank <- apply(pairdist, 1, function(x) malerank[x['d1'], x['d2']])
  pairdist$d2rank <- apply(pairdist, 1, function(x) fmlerank[x['d1'], x['d2']])
  pairdist$distance <- pairdist$d1rank + pairdist$d2rank
  pairdist$correlation <- apply(pairdist, 1, function(x) corsample[x['d1'], x['d2']])
  
  
  ### Determine spurious match due to being left out (threshold is 10% of n)
  #if (sum(1:100 == d1match & pairdist$distance >= 4) > 0){
  spumatch <- which(1:total_samples == d1match & pairdist$distance >= max(2, total_samples/10))
  if (length(spumatch) > 0){
    nonmatch <- c(nonmatch, spumatch)
    cat('  Spurious match:', paste(spumatch), '\n')
  }
  cat('\n Total Number of Mismatched Samples =', length(nonmatch), '\n\n')
  
  
  traincli <- sjcli[, c('sample', clinical_attributes)]
  
  ######### Prediction probability #########
  cli_data_use <- cli_data %>% filter(sample %in% traincli$sample)
  cli_data_use <- cli_data_use[match(traincli$sample, cli_data_use$sample),]

  
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  for(i in 1:length(clinical_attributes)){
    #cli_attr <- clinical_attributes[i]
    cat("clinical attributes: ", clinical_attributes[i], "\n")
    traincli[,clinical_attributes[i]] <- as.factor(traincli[,clinical_attributes[i]])
    traincli[,cli_attr_prob_names_true[i]] <- apply(traincli, 1, function(x) if (x[clinical_attributes[i]] == levels(traincli[,clinical_attributes[i]])[1]) 0 else 1)
  }
  
  
  cat('\n')
  ### Merge probability from teams
  # prediction probability
  cli_attr_prob_names_pred_d1 <- paste("d1",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_d2 <- paste("d2",clinical_attributes,"_prob",sep = "")
  
  cli_attr_prob_names_pred_combine <- paste("pred_",clinical_attributes,sep = "")
  
  for(i in 1:length(clinical_attributes)){
    cat("clinical attributes: ",clinical_attributes[i],"\n")
    # Data_1
    traincli[,cli_attr_prob_names_pred_d1[i]] <- ( sjcli[, cli_attr_prob_names_pred_d1[i]] + sention_d2[, paste(clinical_attributes[i],"_prob",sep = "")] ) / 2.0
    
    # Data_2
    traincli[,cli_attr_prob_names_pred_d2[i]] <- ( sjcli[, cli_attr_prob_names_pred_d2[i]] + sention_d1[, paste(clinical_attributes[i],"_prob",sep = "")] ) / 2.0
    
    # combine Data_1 and Data_2
    traincli[,cli_attr_prob_names_pred_combine[i]] <- ( traincli[,cli_attr_prob_names_pred_d1[i]] + traincli[,cli_attr_prob_names_pred_d2[i]]) / 2.0
  }
  
  
  
  ########### Label Correction
  final_tab <- data.frame(sample=cli_data_use$sample, Clinical=1:total_samples, Data1=1:total_samples, Data2=1:total_samples)
  final_tab <- correctClinicalSwapping(traincli, final_tab, nonmatch, clinical_attributes)
  final_tab <- correctOmicsSwapping(traincli, final_tab, nonmatch, d1match, pairdist, clinical_attributes)
  d1swap    <- which(final_tab$Data1 != 1:total_samples)
  d2swap    <- which(final_tab$Data2 != 1:total_samples)
  swapped   <- c(d1swap, d2swap)
  dup_shift <- setdiff(nonmatch, swapped)
  final_tab <- correctOmicsShifting(traincli, final_tab, dup_shift, corsample, d1match, d2match, swapped, pairdist, clinical_attributes)
  
  
  out_file <- paste(out_dir, "/", prefix, "_final_result.tsv", sep = "")
  write.table(final_tab, out_file, col.names=TRUE, row.names=FALSE, sep='\t')
}




