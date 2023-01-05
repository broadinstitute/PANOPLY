format_input_data = function(pro_file,rna_file,cli_file,out_dir="./"){
    ## load protein and RNA expression data
    proteome <- read.delim(pro_file, stringsAsFactors = FALSE, check.names = FALSE)
    rnaseq   <- read.delim(rna_file, stringsAsFactors = FALSE, check.names = FALSE)
    sam <- read.delim(cli_file,stringsAsFactors = FALSE,check.names = FALSE)
    
    ## use samples present in all omics datasets
    sam_pro_rna <- intersect(names(proteome),names(rnaseq))
    sam_all <- intersect(sam_pro_rna,sam$sample)
    sam_all <- sort(sam_all)
    cat("Remove samples in file ", pro_file,": ", length(setdiff(names(proteome),sam_all)),"\n")
    proteome <- proteome[,sam_all]
    cat("Remove samples in file ", rna_file,": ", length(setdiff(names(rnaseq),sam_all)),"\n")
    rnaseq <- rnaseq[,sam_all]
    sam <- sam[sam$sample %in% sam_all,]
    cat("Use samples:",length(sam_all),"\n")
    
    proteome <- proteome[,sam$sample]
    rnaseq <- rnaseq[,sam$sample]
    
    ## remove rows with all NA
    rows_all_na <- apply(proteome, 1, function(x){all(is.na(x))})
    if(sum(rows_all_na) > 0){
        cat("Remove rows with all NA in file: ",pro_file,", ",sum(rows_all_na),"\n")
        proteome <- proteome[!rows_all_na,]
    }
    
    rows_all_na <- apply(rnaseq, 1, function(x){all(is.na(x))})
    if(sum(rows_all_na) > 0){
        cat("Remove rows with all NA in file: ",rna_file,", ",sum(rows_all_na),"\n")
        rnaseq <- rnaseq[!rows_all_na,]
    }

	cat("Features in file ",pro_file,": ",nrow(proteome),"\n")
	cat("Features in file ",rna_file,": ",nrow(rnaseq),"\n")
	cat("Data range in file ",pro_file,": ",paste(range(proteome,na.rm = TRUE),collapse = " to "),"\n")
	cat("Data range in file ",rna_file,": ",paste(range(rnaseq,na.rm = TRUE),collapse = " to "),"\n")
    
    ## output
    write.table(x = proteome,file = paste(out_dir,"/",basename(pro_file),sep = ""),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
    write.table(x = rnaseq,  file = paste(out_dir,"/",basename(rna_file),sep = ""),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
    write.table(x = sam,     file = paste(out_dir,"/",basename(cli_file),sep = ""),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
}


