library(missForest)
library(biomaRt)
library(glmnet)
library(caret)
library(parallel)
library(doParallel)

set.seed(2020)

run_2b <- function(d1_file, d2_file, anno_file, gene_file, out_dir="./",
                   clinical_attributes=NA){
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ########### Data Retrieval
  cat('Correcting dataset:', d1_file, " and ", d2_file , '... \n')
  
  ## load data
  data1 <- read.delim(d1_file, stringsAsFactors = FALSE, check.names = FALSE)
  data2 <- read.delim(d2_file, stringsAsFactors = FALSE, check.names = FALSE)
  clinic  <- getClinical(anno_file)
  
  ## check for common sample labels
  d1_samples <- colnames(data1)
  d2_samples <- colnames(data2)
  clinic_samples <- clinic[, 'sample']
  
  total_samples  <- union(clinic_samples, union(d1_samples, d2_samples))
  common_samples <- intersect(clinic_samples, intersect(d1_samples, d2_samples))
  sample_n <- length(common_samples)
  if (length(common_samples) < length(total_samples)) {
    cat('Warnings: Found samples with inconsistent label between input data... removed \n')
  }
  data1 <- data1[, common_samples]
  data2 <- data2[, common_samples]
  
  
  out_files <- preprocess(data1, data2, gene_file, out_dir = out_dir)

  data1 <- getOmicsData(out_files[1], out_files[3])
  d1_atsm <- data1$omics_atsm[common_samples, ]
  d1_sex  <- data1$omics_sex[common_samples, ]
  
  data2 <- getOmicsData(out_files[2], out_files[3])
  d2_atsm <- data2$omics_atsm[common_samples, ]
  d2_sex  <- data2$omics_sex[common_samples, ]
  
  
  
  
  
  ########### Screening RNA vs PRO mismatched samples with correlation
  #### first round of getting correlated genes
  cat('Screening First Round: ')
  corsample <- computeCorr(d1_atsm, d2_atsm)
  pairdist <- getMatching(corsample)
  nonmatch <- which(1:sample_n != pairdist$d2)
  cat('  Found', length(nonmatch), 'mismatched samples \n')
  
  #### second round of getting correlated genes (after removing nonmatch from first round)
  cat('Screening Second Round: ')
  corsample <- computeCorr(d1_atsm, d2_atsm, nonmatch=nonmatch)
  pairdist <- getMatching(corsample)
  
  d1match <- pairdist$d2
  d2match <- pairdist$d1[order(pairdist$d2)]
  nonmatch <- which(1:sample_n != d1match)
  cat('  Found', length(nonmatch), 'mismatched samples \n')

  
  # output correlation file
  cor_file <- paste0(out_dir,"/sample_correlation.csv")
  write.table(corsample, cor_file , col.names=TRUE, row.names=TRUE, sep=',')
  
  match_file <- paste0(out_dir, "/pairwise_matching.tsv")
  write.table(pairdist, match_file , row.names=FALSE, col.names=TRUE, sep='\t')
  
  png(paste0(out_dir, '/sample_correlation.png'), width = 600, height = 600)
  heatmap(corsample, Rowv=NA, Colv = NA, scale='none')
  dev.off()
  
  #### Determine spurious match due to being left out (threshold = argmax(2, n*0.1))
  spumatch <- which(1:sample_n == d1match & pairdist$distance >= max(2, sample_n/10))
  if (length(spumatch) > 0){
    nonmatch <- c(nonmatch, spumatch)
    cat('  Spurious match:', paste(spumatch), '\n')
  }
  cat('\n Total Number of mismatched samples =', length(nonmatch), '\n\n')
  
  
  
  ########### Attribute prediction and clinical samples swapping detection
  traincli <- clinic[, c('sample', clinical_attributes)]
  rownames(traincli) <- traincli$sample
  traincli <- traincli[common_samples, ]
  rownames(traincli) <- 1:sample_n
  
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  for(i in 1:length(clinical_attributes)){
    #cli_attr <- clinical_attributes[i]
    cat("clinical attributes: ",clinical_attributes[i],"\n")
    traincli[,clinical_attributes[i]] <- as.factor(traincli[,clinical_attributes[i]])
    traincli[,cli_attr_prob_names_true[i]] <- apply(traincli, 1, function(x) if (x[clinical_attributes[i]] == levels(traincli[,clinical_attributes[i]])[1]) 0 else 1)
  }
  
  #### First round of prediction using cross validation (after removing d1/d2 mismatch samples)
  
  traincli <- predictCV(traincli, nonmatch, d1_sex, d1_atsm, d2_sex, d2_atsm, clinical_attributes)

  cat('\n')
  #### First round flagging potential Clinical Swapping
  #traincli$pred_gender <- (traincli$rgender_prob + traincli$pgender_prob) / 2
  #traincli$pred_msi    <- (traincli$rmsi_prob + traincli$pmsi_prob) / 2
  
  # prediction probability
  cli_attr_prob_names_pred_d1 <- paste("d1",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_d2 <- paste("d2",clinical_attributes,"_prob",sep = "")
  
  cli_attr_prob_names_pred_combine <- paste("pred_",clinical_attributes,sep = "")
  for(i in 1:length(clinical_attributes)){
    cat("clinical attributes: ",clinical_attributes[i],"\n")
    traincli[,cli_attr_prob_names_pred_combine[i]] <- (traincli[,cli_attr_prob_names_pred_d1[i]] + traincli[,cli_attr_prob_names_pred_d2[i]]) / 2.0
  }

  
  clinic_swap <- flagClinicalSwap(traincli, nonmatch, clinical_attributes)
  
  
  #### Second round of prediction (after removing RNA/PRO mismatched samples & clinical swapping cases)
  traincli <- predictLR(traincli, c(nonmatch, clinic_swap), d1_sex, d1_atsm, d2_sex, d2_atsm, clinical_attributes)
  for(i in 1:length(clinical_attributes)){
    cat("clinical attributes: ",clinical_attributes[i],"\n")
    traincli[,cli_attr_prob_names_pred_combine[i]] <- (traincli[,cli_attr_prob_names_pred_d1[i]] + traincli[,cli_attr_prob_names_pred_d2[i]]) / 2.0
  }
  
  write.table(traincli, file = paste(out_dir,"/clinical_attributes_pred.tsv",sep = ""),col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
  
  

  final_tab <- data.frame(sample=rownames(corsample), Clinical=1:sample_n, Data1=1:sample_n, Data2=1:sample_n)
  
  cat('\n')
  ########### Determine Clinical Swapping (Second round to detect and correct the labels)
  final_tab   <- correctClinicalSwapping(traincli, final_tab, nonmatch, clinical_attributes)
  clinic_swap <- which(final_tab$Clinical != 1:sample_n)
  
  cat('\n')
  ########## Determine omics swapping
  final_tab <- correctOmicsSwapping(traincli, final_tab, nonmatch, d1match, pairdist, clinical_attributes)
  rnaswap <- which(final_tab$Data1 != 1:sample_n)
  proswap <- which(final_tab$Data2 != 1:sample_n)
  swapped <- c(rnaswap, proswap)
  
  
  cat('\n')
  ########## Determine omics duplication and shifting
  dup_shift <- setdiff(nonmatch, swapped)
  final_tab <- correctOmicsShifting(traincli, final_tab, dup_shift, corsample, d1match, d2match, swapped, pairdist, clinical_attributes)
  rnashift  <- which(final_tab$Data1 != 1:sample_n)
  rnashift  <- setdiff(rnashift, rnaswap)
  proshift  <- which(final_tab$Data2 != 1:sample_n)
  proshift  <- setdiff(proshift, proswap)
  
  
  ########## Output error statistics and corrected label
  errors <- c(length(clinic_swap), length(proswap), length(proshift), length(rnaswap), length(rnashift))
  names(errors) <- c('cli_swap', 'd2_swap', 'd2_shift', 'd1_swap', 'd1_shift')
  as.data.frame(errors) -> errors
  print(errors)

  cat('\n')
  final_tab_file <- paste0(out_dir, "/final.csv")
  write.table(final_tab, final_tab_file , col.names=TRUE, row.names=FALSE, sep=',')
  error_file <- paste0(out_dir, "/error.tsv")
  write.table(errors, error_file, row.names=T, col.names=T, sep='\t')
}




## Preprocessing: annotate gene with chromosome information
prpr_annotate <- function(geneSymbol, out_gene_file){
  cat('Annotating genes with chromosomes...\n')
  
  gene_info <- tryCatch({
    cat("Try www.ensembl.org ...\n")
    mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "www.ensembl.org")     #ensemblRedirect = FALSE)
    gene_info_tmp <- getBM(attributes=c('hgnc_symbol', 'chromosome_name'),
                       filters='hgnc_symbol',
                       values=geneSymbol,
                       mart=mart)
    return(gene_info_tmp)
  },error=function(e){
    cat("Try http://uswest.ensembl.org/ ...\n")
    mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "http://uswest.ensembl.org/")     #ensemblRedirect = FALSE)
    gene_info_tmp <- getBM(attributes=c('hgnc_symbol', 'chromosome_name'),
                       filters='hgnc_symbol',
                       values=geneSymbol,
                       mart=mart)
    return(gene_info_tmp)
  },
  warning=function(cond){
    print(cond)
    message("Please see the file: kfoldCrossValidation.rda for related data!")
    save(e,ann_use,net_data,r,kfold,ranNum,file="kfoldCrossValidation.rda")
    return(NULL)
  })
  
  cat('Annotation is complete, save annotation into file ', out_gene_file, '\n\n')
  write.table(gene_info, out_gene_file, sep='\t', col.names=T, row.names=F)
}


## preprocessing: omics
prpr_omics <- function(omicsdata, sexgenes) {
  cat('Preprocessing omics data... \n')
  
  ## partition omics data into sex genes and autosomal genes
  sex_omic  <- omicsdata[intersect(rownames(omicsdata), sexgenes), ]
  auto_omic <- omicsdata[setdiff(rownames(omicsdata), sexgenes), ]
  
  # Remove any sex gene which expression values is NA in all samples
  pmiss <- apply(sex_omic, 1, function(x) sum(is.na(x)))
  pmiss <- which(pmiss == ncol(sex_omic))
  if (length(pmiss) > 0) {
    cat('  ', length(pmiss), 'sex gene(s) has NA in ALL samples - removed... \n')
    sex_omic <- sex_omic[-pmiss, ]
  }
  # Replace missing value as 0 in sex genes
  sex_omic <- t(sex_omic)
  sex_omic[is.na(sex_omic)] <- 0
  
  # Remove any automosal gene with > 50% missing values
  pmiss <- apply(auto_omic, 1, function(x) sum(is.na(x)))
  pmiss <- which(pmiss > (ncol(auto_omic)/2))
  if (length(pmiss) > 0) {
    cat('  ', length(pmiss), 'of autosomal gene(s) has NA in > 50% samples - removed... \n')
    auto_omic <- auto_omic[-pmiss, ]
  }
  # Impute missing value in autosomal genes if >= 30% of genes has missing values; else remove
  missingrow <- nrow(auto_omic) - nrow(na.omit(auto_omic))
  missingPct <- missingrow / nrow(auto_omic)
  if (missingPct >= 0.30) {
    cat('  ', sprintf('%d(%.2f%%) of autosomal gene(s) has NA missing values - impute missing values... \n', missingrow, missingPct))
    auto_omic <- t(auto_omic)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    auto_imp <- missForest(auto_omic, parallelize="variables")
    stopCluster(cl)
    auto_omic <- auto_imp$ximp
  } else {
    cat('  ', sprintf('%d(%.2f%%) of autosomal gene(s) has NA missing values - removed... \n', missingrow, missingPct))
    auto_omic <- na.omit(auto_omic)
    auto_omic <- t(auto_omic)
  }
  
  omicsdata <- cbind(sex_omic, auto_omic)
  return(omicsdata)
}


#### Preprocessing Module
preprocess <- function(data1, data2, gene_file, out_dir="./"){
  
  ## annotate genes with chromosomes
  geneSymbol    <- union(rownames(data1), rownames(data2))
  out_gene_file <- paste0(out_dir, "/genes.tsv")
  gene_data     <- read.delim(gene_file, stringsAsFactors = FALSE)
  write.table(gene_data, file = out_gene_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  #prpr_annotate(geneSymbol, out_gene_file)
  
  if(sum(grepl(pattern="^ENSG", x=row.names(data1)))/nrow(data1) > 0.5){
    sexgenes <- getSexGenes(out_gene_file, id_type = "ensg")
  }else{
    sexgenes <- getSexGenes(out_gene_file, id_type = "hgnc_symbol")
  }
  
  ## preprocessing data1
  data1 <- prpr_omics(data1, sexgenes)
  out_data1  <- paste0(out_dir, "/cleaned_data1.tsv")
  write.table(data1, out_data1, row.names=TRUE, col.names=TRUE, sep='\t')
  
  ## preprocessing data2
  data2   <- prpr_omics(data2, sexgenes)
  out_data2  <- paste0(out_dir, "/cleaned_data2.tsv")
  write.table(data2, out_data2, row.names=TRUE, col.names=TRUE, sep='\t')
  
  return(c(out_data1, out_data2, out_gene_file))
}




## Load: gene chromosome annotation
getSexGenes <- function(gene_file, id_type="hgnc_symbol"){
  gene_info <- read.delim(gene_file, stringsAsFactors = FALSE)
  sexgenes <- which(gene_info$chromosome_name %in% c('X', 'Y'))
  sexgenes <- unique(gene_info[, id_type][sexgenes])
  return(sexgenes)
}


## Load sample annotation data: test_cli.tsv
getClinical <- function(anno_file){
  clinic  <- read.delim(anno_file, stringsAsFactors = FALSE,check.names = FALSE)
  return(clinic)
}


## Load omics data
getOmicsData <- function(d_file, gene_file){
  omicsdata   <- read.delim(d_file, stringsAsFactors = FALSE)
  
  # each column is a gene
  if(sum(grepl(pattern="^ENSG",x=colnames(omicsdata)))/ncol(omicsdata) > 0.5){
    sexgenes <- getSexGenes(gene_file,id_type = "ensg")
  }else{
    sexgenes <- getSexGenes(gene_file,id_type = "hgnc_symbol")
  }
  
  #sexgenes <- getSexGenes(gene_file)
  
  omics_atsm <- omicsdata[, setdiff(colnames(omicsdata), sexgenes)]
  omics_sex  <- omicsdata[, intersect(colnames(omicsdata), sexgenes)]
  
  omicsdata <- list()
  omicsdata$omics_atsm <- scale(omics_atsm)
  omicsdata$omics_sex  <- as.matrix(omics_sex)
  return(omicsdata)
}





## Screen: compute correlation matrix between d1&d2 samples
computeCorr <- function(d1_atsm, d2_atsm, nonmatch = c()) {
  intergene <- intersect(colnames(d1_atsm), colnames(d2_atsm))
  if (length(nonmatch) == 0) {
    corgene <- extractCorGene(d1_atsm[, intergene], d2_atsm[, intergene])
  } else {
    corgene <- extractCorGene(d1_atsm[-nonmatch, intergene], d2_atsm[-nonmatch, intergene])
  }
  corsample <- cor(t(d1_atsm[, corgene]), t(d2_atsm[, corgene]))
  return(corsample)
}


## Screen: extract correlated genes between d1&d2 omics data
extractCorGene <- function(d1matrix, d2matrix){
  cormatrix <- cor(as.matrix(d1matrix), as.matrix(d2matrix))
  correlate <- diag(cormatrix)
  corgene <- names(which(correlate > 0.5))
  cat(sprintf('%d genes are highly correlated \n', length(corgene)))
  return(corgene)
}


## Screen: convert correlation matrix into probability matrix
computeRankdist <- function(corsample){
  dist_d1_d2 <- corsample
  dist_d2_d1 <- corsample
  for (i in 1:nrow(corsample)){
    dist_d1_d2[i,] <- exp(scale(dist_d1_d2[i,])) / sum(exp(scale(dist_d1_d2[i,])))
    dist_d2_d1[,i] <- exp(scale(dist_d2_d1[,i])) / sum(exp(scale(dist_d2_d1[,i])))
  }
  rankdist <- dist_d1_d2 * dist_d2_d1
  
  return(rankdist)
}


## Screen: Perform Stable Matching between omics samples and return the matching with scores
getMatching <- function(cormatrix){
  row_labels <- rownames(cormatrix)
  col_labels <- colnames(cormatrix)
  rownames(cormatrix) <- paste0('d1_', 1:nrow(cormatrix))
  colnames(cormatrix) <- paste0('d2_', 1:ncol(cormatrix))
  
  probmatrix <- computeRankdist(cormatrix)
  malerank <- rankbyRow(probmatrix)
  fmlerank <- rankbyCol(probmatrix)
  
  matcher <- stableMarriage(malerank, fmlerank)
  d1match <- as.numeric(sub('d2_', '', matcher$malematch))
  d2match <- as.numeric(sub('d1_', '', matcher$fmlematch))
  
  pairdist <- data.frame(d1=1:length(d1match), d1_label=row_labels, d2=d1match, d2_label=col_labels[d1match])
  pairdist$d1rank <- apply(pairdist, 1, function(x) malerank[as.numeric(x['d1']), as.numeric(x['d2'])])
  pairdist$d2rank <- apply(pairdist, 1, function(x) fmlerank[as.numeric(x['d1']), as.numeric(x['d2'])])
  pairdist$distance <- pairdist$d1rank + pairdist$d2rank
  pairdist$correlation <- apply(pairdist, 1, function(x) cormatrix[as.numeric(x['d1']), as.numeric(x['d2'])])
  
  return(pairdist)
}


## StableMatching: get preferential rank by row
rankbyRow <- function(rankdist){
  malerank <- t(apply(rankdist, 1, function(x) rank(-x, ties.method='first')))
  colnames(malerank) <- colnames(rankdist)
  return(malerank)
}


## StableMatching: get preferential rank by column
rankbyCol <- function(rankdist) {
  fmlerank <- apply(rankdist, 2, function(x) rank(-x, ties.method='first'))
  rownames(fmlerank) <- rownames(rankdist)
  return(fmlerank)
}


## StableMatching: get matching pairs of samples with preferential ranks
stableMarriage <- function(malerank, fmlerank){
  matcher <- list()
  
  malematch <- rep(NA, nrow(malerank))
  names(malematch) <- rownames(malerank)
  fmlematch <- rep(NA, ncol(fmlerank))
  names(fmlematch) <- colnames(fmlerank)
  
  singlemales <- names(malematch)
  while (length(singlemales) != 0){
    for (ppsmale in singlemales){
      propose <- 1
      single <- TRUE
      while (single == TRUE){
        ppsfmle <- names(which(malerank[ppsmale,] == propose))
        engaged <- fmlematch[ppsfmle]
        if (is.na(engaged) || fmlerank[engaged, ppsfmle] > fmlerank[ppsmale, ppsfmle]) {
          malematch[ppsmale] <- ppsfmle
          fmlematch[ppsfmle] <- ppsmale
          singlemales <- setdiff(singlemales, ppsmale)
          if (!(is.na(engaged))) singlemales <- c(singlemales, engaged) 
          single <- FALSE
        } else {
          propose <- propose + 1
        }
      }
    }
  }
  
  matcher$malematch <- malematch
  matcher$fmlematch <- fmlematch
  return(matcher)
}




## Train: get weight of training instance, distributed equally by class
getClassWeight <- function(labels){
  weight <- rep( 1, length(labels) )
  labelclass <- unique(labels)
  weight[labels == labelclass[1]] <- length(labels)/ 2 / sum(labels == labelclass[1])
  weight[labels == labelclass[2]] <- length(labels)/ 2 / sum(labels == labelclass[2])
  return(weight)
}


## Train: weighted L1 regularized Logistic Regression, using CV to determine best lambda
trainGLM <- function(msiLabel, inputmtx, alpha){
  weight <- getClassWeight(msiLabel)
  if (sum(weight) != length(msiLabel)){
    cat('Error: sum(classweight) does not equal to length(msiLabel)! \n')
  }   # should be nrow(trainset)
  
  # perform cross validation of elasticnet to determine optimum lambda
  
  cv.glm <- cv.glmnet(as.matrix(inputmtx), msiLabel, family="binomial", weights=weight, alpha=alpha, parallel=TRUE)
  (best_lambda <- cv.glm$lambda.1se)
  fit <- glmnet(as.matrix(inputmtx), msiLabel, family="binomial", weights=weight, alpha=alpha, lambda=best_lambda)
  return(fit)
}


## Train: use trainGLM fn in k-fold
trainGLMcv <- function(msiLabel, inputmtx, alpha, k = 5){
  flds <- createFolds(msiLabel, k = k, list = TRUE, returnTrain = FALSE)
  predoutput <- data.frame(att=msiLabel, att_prob=0)
  
  for (f in 1:k) {
    testidx <- flds[[f]]
    numvar <- 0
    iter   <- 0
    cat('  Training Fold', f, '- Optimizing Model...')
    while (numvar < 4 && iter < 50){
      
      fit1 <- tryCatch({
        fit1 <- trainGLM(msiLabel[-testidx], inputmtx[-testidx, ], 0.3)
        fit1
      },error=function(e){
        save(msiLabel,testidx,inputmtx,testidx,f,file="trainGLM_input_data.rda")
        cat("\n------\n")
        print(e)
        cat("\n------\n")
        stop("Error in trainGLM\n")
        return(NULL)
      })
      
      numvar1 <- sum(coef(fit1) > 0)
      if (numvar1 >= numvar) {
        fit <- fit1
        numvar <- numvar1
      }
      iter <- iter + 1
    }
    cat('Done. \n')
    predoutput$att[testidx] <- predict(fit, as.matrix(inputmtx[testidx, ]), type='class')[, 1]
    predoutput$att_prob[testidx] <- predict(fit, as.matrix(inputmtx[testidx, ]), type='response')[, 1]
  }
  return(predoutput)
}

is_gender = function(x){
  gd <- grepl(pattern = "gender",x = x, ignore.case = TRUE)
  return(gd)
}


## Train: use trainGLMcv fn for each clinical attribute
# clinical_attributes = c("gender")
predictCV <- function(traincli, nonmatch, d1_sex, d1_atsm, d2_sex, d2_atsm, 
                      clinical_attributes){
  
  matchingidx <- setdiff(1:nrow(traincli), nonmatch)
  
  # prediction probability
  cli_attr_prob_names_pred_d1 <- paste("d1",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_d2 <- paste("d2",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_d1 <- paste("d1",clinical_attributes,sep = "")
  cli_attr_prob_names_d2 <- paste("d2",clinical_attributes,sep = "")
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  for(i in 1:length(clinical_attributes)){
    if(is_gender(clinical_attributes[i])){
      ## gender prediction
      cat("Predicting ", clinical_attributes[i] ," from data1... \n")
      save(traincli,clinical_attributes,i,nonmatch,d1_sex,d2_sex,file="trainGLMcv_input_debug.rda")
      predoutput_d1 <- trainGLMcv(traincli[,clinical_attributes[i]][matchingidx], d1_sex[matchingidx, ], 0.3)
      cat("Predicting ", clinical_attributes[i] ," from data2... \n")
      predoutput_d2 <- trainGLMcv(traincli[,clinical_attributes[i]][matchingidx], d2_sex[matchingidx, ], 0.3)
      
    }else{
      ## non-gender attribute prediction
      cat("Predicting ", clinical_attributes[i] ," from data1... \n")
	    save(traincli,clinical_attributes,i,nonmatch,d1_atsm,d2_atsm,file="trainGLMcv_input_debug.rda")
      predoutput_d1 <- trainGLMcv(traincli[,clinical_attributes[i]][matchingidx], d1_atsm[matchingidx, ], 0.3)
      cat("Predicting ", clinical_attributes[i] ," from data2... \n")
      predoutput_d2 <- trainGLMcv(traincli[,clinical_attributes[i]][matchingidx], d2_atsm[matchingidx, ], 0.3)
    }
    
    ## data1
    traincli[,cli_attr_prob_names_d1[i]] <- traincli[,clinical_attributes[i]]
    traincli[,cli_attr_prob_names_pred_d1[i]] <- traincli[,cli_attr_prob_names_true[i]]
    traincli[,cli_attr_prob_names_d1[i]][matchingidx] <- predoutput_d1$att
    traincli[,cli_attr_prob_names_pred_d1[i]][matchingidx] <- predoutput_d1$att_prob
    
    ## data2
    traincli[,cli_attr_prob_names_d2[i]] <- traincli[,clinical_attributes[i]]
    traincli[,cli_attr_prob_names_pred_d2[i]] <- traincli[,cli_attr_prob_names_true[i]]
    traincli[,cli_attr_prob_names_d2[i]][matchingidx] <- predoutput_d2$att
    traincli[,cli_attr_prob_names_pred_d2[i]][matchingidx] <- predoutput_d2$att_prob
    
  }
  
  stopCluster(cl)
  
  return(traincli)
}


find_gender_label = function(clinical_attributes = c("gender")){
  
  gender_label <- NA
  for(i in clinical_attributes){
    if(is_gender(i)){
      gender_label <- i
    }
  }
  return(gender_label)
}


## Correction: Detect potential clinical swapping
flagClinicalSwap <- function(traincli, nonmatch, clinical_attributes) {
  
  gender_label <- find_gender_label(clinical_attributes)
  if(!is.na(gender_label)){
    gender_prob <- paste(gender_label,"_prob",sep = "")
    pred_gender <- paste("pred_",gender_label,sep = "")
    #high_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.7)
    high_suspect <- which(abs(traincli[,gender_prob] - traincli[,pred_gender]) > 0.7)
    (high_suspect <- setdiff(high_suspect, nonmatch))
    cli_suspect <- which(abs(traincli[,gender_prob] - traincli[,pred_gender]) > 0.5)
    (cli_suspect <- setdiff(cli_suspect, nonmatch))
  }else{
    cli_suspect <- c()
    high_suspect <- c() # SV ADDED THIS LINE IN!
  }
  
  # true probability
  cli_attr_prob_names_true <- paste0(clinical_attributes, "_prob")
  cli_attr_prob_names_pred_combine <- paste0("pred_", clinical_attributes)
  
  
  clinic_swap <- c()
  if (length(cli_suspect) <= 1){
    cat('No Clinical Swapping Cases Found! \n')
  } else  {
    #subset <- traincli[cli_suspect, c('gender_prob', 'msi_prob', 'pred_gender', 'pred_msi')]
    subset <- traincli[cli_suspect, c(cli_attr_prob_names_true, cli_attr_prob_names_pred_combine)]
    
    
    clidist <- matrix(nrow=length(cli_suspect), ncol=length(cli_suspect))
    rownames(clidist) <- cli_suspect
    colnames(clidist) <- cli_suspect
    for (male in as.character(cli_suspect)) {
      for (fmle in as.character(cli_suspect)) {
        #clidist[male, fmle] <- (1 - abs(subset[male, 'gender_prob'] - subset[fmle, 'pred_gender'])) + (1 - abs(subset[male, 'msi_prob'] - subset[fmle, 'pred_msi']))
        clidist[male, fmle] <- 0
        for(i in 1:length(clinical_attributes)){
          clidist[male, fmle] <- clidist[male, fmle] + (1 - abs(subset[male, cli_attr_prob_names_true[i]] - subset[fmle, cli_attr_prob_names_pred_combine[i]]))
        }
      }
    }
    
    climale <- rankbyRow(clidist)
    clifmle <- t(climale)
    
    clinic_match <- stableMarriage(climale, clifmle)
    for (eachmatch in names(clinic_match$malematch)){
      match1 <- as.numeric(eachmatch)
      match2 <- as.numeric(clinic_match$malematch[eachmatch])
      if (match1 %in% clinic_swap)    next
      if (match1 != match2) {
        clinic_swap <- c(clinic_swap, match1, match2)
      }
    }
  }
  clinic_swap <- union(clinic_swap, high_suspect)
  return(clinic_swap)
}


## Train: use trainGLM fn for each clinical attribute
# clinical_attributes = c("gender")
predictLR <- function(traincli, nonmatch, d1_sex, d1_atsm, d2_sex, d2_atsm,
                      clinical_attributes) {
  
  matchingidx <- setdiff(1:nrow(traincli), nonmatch)
  
  # prediction probability
  cli_attr_prob_names_pred_d1 <- paste("d1",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_d2 <- paste("d2",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_d1 <- paste("d1",clinical_attributes,sep = "")
  cli_attr_prob_names_d2 <- paste("d2",clinical_attributes,sep = "")
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  for(i in 1:length(clinical_attributes)){
    if(is_gender(clinical_attributes[i])){
      numvar <- 0
      iter <- 0
      cat('Training data1 model for ',clinical_attributes[i],"\n",sep = "")
      while (numvar < 4 && iter < 50){
        fit <- trainGLM(traincli[,clinical_attributes[i]][matchingidx], d1_sex[matchingidx, ], 0.3)
        numvar <- sum(coef(fit) > 0)
        iter <- iter + 1
      }
      cat('Done - Regularized model with', numvar, 'variables. \n')
      traincli[,cli_attr_prob_names_d1[i]] <- predict(fit, d1_sex, type='class')[, 1]
      traincli[,cli_attr_prob_names_pred_d1[i]] <- predict(fit, d1_sex, type='response')[, 1]
      
      numvar <- 0
      iter   <- 0
      cat('Training data2 model for ',clinical_attributes[i],"\n",sep = "")
      while (numvar < 4 && iter < 50){
        fit1 <- trainGLM(traincli[,clinical_attributes[i]][matchingidx], d2_sex[matchingidx, ], 0.3)
        numvar1 <- sum(coef(fit1) > 0)
        if (numvar1 > numvar) {
          fit <- fit1
          numvar <- numvar1
        }
        iter <- iter + 1
      }
      cat('Done - Regularized model with', numvar, 'variables. \n')
      traincli[,cli_attr_prob_names_d2[i]] <- predict(fit, d2_sex, type='class')[, 1]
      traincli[,cli_attr_prob_names_pred_d2[i]] <- predict(fit, d2_sex, type='response')[, 1]
      
      
    } else {
      ## non-gender attribute prediction
      numvar <- 0
      iter <- 0
      cat('Training data1 model for ',clinical_attributes[i],"\n",sep = "")
      while (numvar < 4){
        fit <- trainGLM(traincli[,clinical_attributes[i]][matchingidx], d1_atsm[matchingidx, ], 0.3)
        numvar <- sum(coef(fit) > 0)
        iter <- iter + 1
      }
      cat('Done - Regularized model with', numvar, 'variables. \n')
      traincli[,cli_attr_prob_names_d1[i]] <- predict(fit, d1_atsm, type='class')[, 1]
      traincli[,cli_attr_prob_names_pred_d1[i]] <- predict(fit, d1_atsm, type='response')[, 1]
      
      numvar <- 0
      iter   <- 0
      cat('Training data2 model for ',clinical_attributes[i],"\n",sep = "")
      while (numvar < 4 && iter < 50){
        fit1 <- trainGLM(traincli[,clinical_attributes[i]][matchingidx], d2_atsm[matchingidx, ], 0.3)
        numvar1 <- sum(coef(fit1) > 0)
        if (numvar1 > numvar) {
          fit <- fit1
          numvar <- numvar1
        }
        iter <- iter + 1
      }
      cat('Done - Regularized model with', numvar, 'variables. \n')
      traincli[,cli_attr_prob_names_d2[i]] <- predict(fit, d2_atsm, type='class')[, 1]
      traincli[,cli_attr_prob_names_pred_d2[i]] <- predict(fit, d2_atsm, type='response')[, 1]
      
    }
    
  }
  
  stopCluster(cl)
  
  
  return(traincli)
}


## Correction: Detect clinical swapping and correct label
correctClinicalSwapping <- function(traincli, final_tab, nonmatch, clinical_attributes) {

  gender_label <- find_gender_label(clinical_attributes)
  if(!is.na(gender_label)){
    gender_prob <- paste(gender_label,"_prob",sep = "")
    pred_gender <- paste("pred_",gender_label,sep = "")
    #high_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.7)
    high_suspect <- which(abs(traincli[,gender_prob] - traincli[,pred_gender]) > 0.7)
    (high_suspect <- setdiff(high_suspect, nonmatch))
    if (length(high_suspect) > 0) {
      cat('Highly suspected Clinical swap case:', paste(high_suspect), '\n')
      final_tab$Clinical[high_suspect] <- -1
    }
    cli_suspect <- which(abs(traincli[,gender_prob] - traincli[,pred_gender]) > 0.55)
    cli_suspect <- setdiff(cli_suspect, nonmatch)
  }else{
    cli_suspect <- c()
  }
  
  
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_combine <- paste("pred_",clinical_attributes,sep = "")
  
  clinic_swap <- c()
  if (length(cli_suspect) <= 1){
    cat('No Clinical Swapping Cases Found! \n')
  } else  {
    #subset <- traincli[cli_suspect, c('gender_prob', 'msi_prob', 'pred_gender', 'pred_msi')]
    subset <- traincli[cli_suspect, c(cli_attr_prob_names_true,cli_attr_prob_names_pred_combine)]
    
    clidist <- matrix(nrow=length(cli_suspect), ncol=length(cli_suspect))
    rownames(clidist) <- cli_suspect
    colnames(clidist) <- cli_suspect
    for (male in as.character(cli_suspect)) {
      for (fmle in as.character(cli_suspect)) {
        #clidist[male, fmle] <- (1 - abs(subset[male, 'gender_prob'] - subset[fmle, 'pred_gender'])) + (1 - abs(subset[male, 'msi_prob'] - subset[fmle, 'pred_msi']))
        clidist[male, fmle] <- 0
        for(i in 1:length(clinical_attributes)){
          clidist[male, fmle] <- clidist[male, fmle] + (1 - abs(subset[male, cli_attr_prob_names_true[i]] - subset[fmle, cli_attr_prob_names_pred_combine[i]]))
        }
      }
    }
    
    climale <- rankbyRow(clidist)
    clifmle <- t(climale)
    
    clinic_match <- stableMarriage(climale, clifmle)
    for (eachmatch in names(clinic_match$malematch)){
      match1 <- as.numeric(eachmatch)
      match2 <- as.numeric(clinic_match$malematch[eachmatch])
      if (match1 %in% clinic_swap)    next
      if (match1 != match2) {
        final_tab[match1, 'Clinical'] <- match2
        final_tab[match2, 'Clinical'] <- match1
        clinic_swap <- c(clinic_swap, match1, match2)
        cat(sprintf('Clinical Label Swapping: %d <--> %d \n', match1, match2))
      }
    }
    print(traincli[clinic_swap, ])
  }

  return(final_tab)
}


## Correct: Detect omics swapping and correct label
correctOmicsSwapping <- function(traincli, final_tab, nonmatch, d1match, pairdist,
                                 clinical_attributes) {
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  # prediction probability
  cli_attr_prob_names_pred_d1 <- paste("d1",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_d2 <- paste("d2",clinical_attributes,"_prob",sep = "")
  
  swapped <- c()
  for (r in nonmatch){
    if (r %in% swapped)    next
    p <- d1match[r]
    if (d1match[p] == r && d1match[r] != r) {
      if (sum(pairdist[c(r, p), c('d1rank', 'd2rank')]) > max(2, nrow(traincli)*0.1) ){
        cat(sprintf('Spurious Pair: %d <--> %d, distance = %f \n', r, p, sum(pairdist[c(r, p), c('d1rank', 'd2rank')])))
        next
      }
      
      swapped <- c(swapped, r, p)
      #subset <- traincli[c(r,p), c('gender_prob', 'rgender_prob', 'pgender_prob', 'msi_prob', 'rmsi_prob', 'pmsi_prob')]
      #pro_swap <- abs(subset[1,1]-subset[1,2]) + abs(subset[1,1]-subset[2,3]) + abs(subset[2,1]-subset[2,2]) + abs(subset[2,1]-subset[1,3])
      #pro_swap <- pro_swap + abs(subset[1,4]-subset[1,5]) + abs(subset[1,4]-subset[2,6]) + abs(subset[2,4]-subset[2,5]) + abs(subset[2,4]-subset[1,6])
      #rna_swap <- abs(subset[1,1]-subset[2,2]) + abs(subset[1,1]-subset[1,3]) + abs(subset[2,1]-subset[1,2]) + abs(subset[2,1]-subset[2,3])
      #rna_swap <- rna_swap + abs(subset[1,4]-subset[2,5]) + abs(subset[1,4]-subset[1,6]) + abs(subset[2,4]-subset[1,5]) + abs(subset[2,4]-subset[2,6])
      
      pro_swap <- 0
      rna_swap <- 0
      for(i in 1:length(clinical_attributes)){
        subset <- traincli[c(r,p), c(cli_attr_prob_names_true[i], cli_attr_prob_names_pred_d1[i], cli_attr_prob_names_pred_d2[i])]
        pro_swap <- pro_swap + abs(subset[1,1]-subset[1,2]) + abs(subset[1,1]-subset[2,3]) + abs(subset[2,1]-subset[2,2]) + abs(subset[2,1]-subset[1,3])
        rna_swap <- rna_swap + abs(subset[1,1]-subset[2,2]) + abs(subset[1,1]-subset[1,3]) + abs(subset[2,1]-subset[1,2]) + abs(subset[2,1]-subset[2,3])
      }
      
      
      if (pro_swap < rna_swap){
        cat(sprintf('Data2 swap: %d <--> %d (d1_error: %.3f vs d2_error: %.3f) \n', r, p, rna_swap, pro_swap))
        final_tab[r, 'Data2'] <- p
        final_tab[p, 'Data2'] <- r
      } else {
        cat(sprintf('Data1 swap: %d <--> %d (d1_error: %.3f vs d2_error: %.3f) \n', r, p, rna_swap, pro_swap))
        final_tab[r, 'Data1'] <- p
        final_tab[p, 'Data1'] <- r
      }
    }
  }
  return(final_tab)
}


## Correct: Detect omics duplication + shifting and correct label
correctOmicsShifting <- function(traincli, final_tab, dup_shift, cormatrix, 
                                 d1match, d2match, swapped, pairdist,
                                 clinical_attributes) {
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  # prediction probability
  cli_attr_prob_names_pred_d1 <- paste("d1",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_d2 <- paste("d2",clinical_attributes,"_prob",sep = "")
  
  rankdist <- computeRankdist(cormatrix)
  rnashift <- c()
  proshift <- c()
  if (length(dup_shift) == 0) {
    cat('No shifting & duplication suspect! Skip... \n')
  } else {
    cat(length(dup_shift), 'samples suspected for duplication and shifting:', paste(dup_shift), '\n')
    shiftdist <- pairdist[dup_shift,]
    shiftdist <- shiftdist[order(-shiftdist$distance), ]
    print(shiftdist)
    
    lose_starts <- shiftdist$d1[shiftdist$distance > 3]
    lose_ends   <- shiftdist$d2[shiftdist$distance > 3]
    
    
    ### chain identification
    chains <- list()
    i <- 1
    for (start in lose_starts){
      chain <- c(start)
      cnext <- start
      while (!(cnext %in% lose_ends)) {
        cnext <- d2match[cnext]
        chain <- c(chain, cnext)
      }
      
      chain <- c(chain, which.max(rankdist[, chain[length(chain)]]))
      chain <- c(which.max(rankdist[chain[1],]), chain)
      
      chains[[i]] <- chain
      i <- i + 1
    }
    
    if (sum(!(dup_shift %in% unlist(chains))) == 0){
      cat('All suspected samples are found in chain! \n\n')
    } else {
      cat('Warning! Not all suspected samples are found in chain. Circular shifting suspected! \n\n')
    }
    
    
    for (chain in chains){
      cat(paste(chain, collapse = ' --> '), '\n')
      if (chain[1] %in% swapped) {
        cat('Warning: chain head', chain[1], 'is swapped samples!\n')
      }
      if (chain[length(chain)] %in% swapped) {
        cat('Warning: chain tail', chain[length(chain)], 'is swapped samples!\n')
      }
      
      #subset <- traincli[chain, c('gender_prob', 'rgender_prob', 'pgender_prob', 'msi_prob', 'rmsi_prob', 'pmsi_prob')]
      #lenchain <- length(chain)
      #rna_shift <- sum(abs(subset[1:(lenchain-2), 1] - subset[2:(lenchain-1), 2])) + sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 3])) + abs(subset[1,1] - subset[1,2])
      #rna_shift <- rna_shift + sum(abs(subset[1:(lenchain-2), 4] - subset[2:(lenchain-1), 5])) + sum(abs(subset[2:lenchain, 4] - subset[2:lenchain, 6])) + abs(subset[1,4] - subset[1,5])
      #pro_shift <- sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 2])) + sum(abs(subset[3:lenchain, 1] - subset[2:(lenchain-1), 3])) + abs(subset[lenchain,1] - subset[lenchain,3])
      #pro_shift <- pro_shift + sum(abs(subset[2:lenchain, 4] - subset[2:lenchain, 5])) + sum(abs(subset[3:lenchain, 4] - subset[2:(lenchain-1), 6])) + abs(subset[lenchain,4] - subset[lenchain,6])
      
      rna_shift <- 0
      pro_shift <- 0
      lenchain <- length(chain)
      for(i in 1:length(clinical_attributes)){
        subset <- traincli[chain, c(cli_attr_prob_names_true[i], cli_attr_prob_names_pred_d1[i], cli_attr_prob_names_pred_d2[i])]
        rna_shift <- rna_shift + sum(abs(subset[1:(lenchain-2), 1] - subset[2:(lenchain-1), 2])) + sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 3])) + abs(subset[1,1] - subset[1,2])
        pro_shift <- pro_shift + sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 2])) + sum(abs(subset[3:lenchain, 1] - subset[2:(lenchain-1), 3])) + abs(subset[lenchain,1] - subset[lenchain,3])
      }
      ## if rna_shift < pro_shift, means error rate for rna shifting is lower, means it is rna shifting
      
      distfront <- rankdist[chain[2], chain[1]]
      distback  <- rankdist[chain[lenchain], chain[lenchain-1]]
      ## if distback > distfront, means it is likely proteome duplication than RNAseq duplication
      
      mean_true <- c()
      for(i in 1:length(clinical_attributes)){
        subset <- traincli[chain, c(cli_attr_prob_names_true[i], cli_attr_prob_names_pred_d1[i], cli_attr_prob_names_pred_d2[i])]
        mean_true[i] <- mean(subset[, cli_attr_prob_names_true[i]][2:(lenchain-1)]) == 0 || mean(subset[, cli_attr_prob_names_true[i]][2:(lenchain-1)]) == 1
      }
      
      
      ## if all samples in a shifting chain the have same attribute (e.g. consistenly shifting ALL MALE & Low-MSI sample)
      #if ((mean(subset$gender_prob[2:(lenchain-1)]) == 0 || mean(subset$gender_prob[2:(lenchain-1)]) == 1) && (mean(subset$msi_prob[2:(lenchain-1)]) == 0 || mean(subset$msi_prob[2:(lenchain-1)]) == 1)) {
      if(all(mean_true)){
        if (distback > distfront) {
          cat('Same attr, so Data2 shift: ', chain[1], paste(chain[3:lenchain], collapse=' '), chain[length(chain)], sprintf('(d1_error: %.3f vs d2_error: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: d1_%d <-> d2_%d = %.4f \t d1_%d <-> d2_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'Data2'] <- chain[3:length(chain)]
          proshift <- c(proshift, chain[2:(lenchain-1)])
        } else {
          cat('Same attr, so Data1 shift: ', chain[1], paste(chain[1:(length(chain)-2)], collapse = ' '), chain[length(chain)], sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: d1_%d <-> d2_%d = %.4f \t d1_%d <-> d2_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'Data1'] <- chain[1:(length(chain)-2)]
          rnashift <- c(rnashift, chain[2:(lenchain-1)])
        }
        ## if all samples in a shifting chain the have different attribute (e.g. shifting of M/F or Low/high samples)
      } else{
        if (pro_shift < rna_shift) {          # supposedly distfront < distback
          cat('Data2 shift: ', chain[1], paste(chain[3:lenchain], collapse=' '), chain[length(chain)], sprintf('(d1_error: %.3f vs d2_error: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: d1_%d <-> d2_%d = %.4f \t d1_%d <-> d2_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'Data2'] <- chain[3:length(chain)]
          proshift <- c(proshift, chain[2:(lenchain-1)])
          if (distfront > distback) {
            cat('Warning: distance indicates Data1_', chain[1], ' duplication but prediction results indicate Data2 shifting', sprintf('(d1_error: %.3f vs d2_error: %.3f) \n', rna_shift, pro_shift), '\n')
          }
          
        } else if (rna_shift < pro_shift) {   # supposedly distfront > distback
          cat('Data1 shift: ', chain[1], paste(chain[1:(length(chain)-2)], collapse = ' '), chain[length(chain)], sprintf('(d1_error: %.3f vs d2_error: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: d1_%d <-> d2_%d = %.4f \t d1_%d <-> d2_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'Data1'] <- chain[1:(length(chain)-2)]
          rnashift <- c(rnashift, chain[2:(lenchain-1)])
          if (distfront < distback){
            cat('Warning: distance indicates Data2_', chain[lenchain], ' duplication but prediction results indicate Data1 shifting', sprintf('(d1_error: %.3f vs d2_error: %.3f) \n', rna_shift, pro_shift), '\n')
          }
        }
      }
      cat('\n')
    }
  }
  
  return(final_tab)
}


