library(blacksheepr)
library(cmapR)
library(dplyr)
library(tibble)
library(stringr)
library(yaml)

args <- commandArgs(TRUE)

gct_path = as.character(args[1])
yaml_file = as.character(args[2])
#groups_file = as.character(args[3])

yaml_params = read_yaml(yaml_file)
GeneSymbol_column = yaml_params$global_parameters$gene_mapping$gene_id_col
identifiers_file = yaml_params$panoply_blacksheep$identifiers_file
groups_file = yaml_params$panoply_blacksheep$groups_file
fraction_samples_cutoff = yaml_params$panoply_blacksheep$fraction_samples_cutoff
fdrcutoffvalue = yaml_params$panoply_blacksheep$fdr_value

create_values_input = function(gct, GeneSymbol_column, identifiers_file){
  
  # extract data table and replace protein name with gene symbol for downstream aggregation and analysis
  data_values = data_values %>%
    rownames_to_column("rowname")

  genesymbol = gct@rdesc %>%
    rownames_to_column("rowname") %>%
    select(all_of(GeneSymbol_column), rowname) %>%
    left_join(data_values, by = "rowname")
  
  if (identifiers_file != NULL){
    if (identifiers_file == "kinases"){
      identifiers_file = "/prot/proteomics/Projects/PGDAC/src/kinase_list.txt"
    }
    identifiers = read.delim(identifiers_file, header = FALSE, sep="\t")
    identifiers = as.character(identifiers$V1)
    genesymbol = genesymbol %>%
      filter(genesymbol[,GeneSymbol_column] %in% identifiers)
  }
  
  genesymbol[,GeneSymbol_column] = paste(genesymbol[,GeneSymbol_column], row.names(genesymbol), sep = "-")
  
  genesymbol_values = genesymbol %>%
    column_to_rownames(GeneSymbol_column) %>%
    select(-rowname)
  
  return(genesymbol_values)
}

create_annotations_input = function(gct, groups_file){
  
  # extract annotations of interest for outlier analysis groupings
  annotations = gct@cdesc %>%
    select(all_of(colnames(groups_file)))
  
  # binarize annotations
  binary_annotations = make_comparison_columns(annotations)
  
  return(binary_annotations)
}

outlier_count_only = function(data_values){

  csv_name = c("negative", "positive")
  x = c(TRUE,FALSE)
  
  for (i in 1:length(x)){
    
    dir.create(csv_name[i])
    
    # make outliers table
    outlier_table_out = make_outlier_table(data_values,
                                           analyze_negative_outliers = x[i])
    
    # save table indicating outliers
    outlier_table = outlier_table_out$outliertab
    write.csv(outlier_table, file.path(csv_name[i], paste0(csv_name[i], "_outliers_table.csv")))

  return(outlier_table
}

outlier_count_and_analysis = function(groupings, data_values, fraction_samples_cutoff){
  results = vector("list", 2)
  csv_name = c("negative", "positive")
  names(results) = csv_name
  x = c(TRUE,FALSE)
  
  for (i in 1:length(x)){
    
    dir.create(csv_name[i])
    
    # make outliers table
    outlier_table_out = make_outlier_table(data_values,
                                           analyze_negative_outliers = x[i])
    
    # save table indicating outliers
    outlier_table = outlier_table_out$outliertab
    write.csv(outlier_table, file.path(csv_name[i], paste0(csv_name[i], "_outliers_table.csv")))
    
    # tabulate outliers per feature
    count_outliers_out <- count_outliers(groupings, outlier_table, 
                                         aggregate_features = TRUE, feature_delineator = "-")
    
    # save aggregated outlier table - number of outliers per feature (gene symbol)
    write.csv(count_outliers_out$aggoutliertab, file.path(csv_name[i], paste0("aggregated_", csv_name[i], "_outliers_table.csv")))
    
    # save fraction table - fraction of outliers per aggregated feature per sample
    fraction_table = count_outliers_out$fractiontab
    write.csv(fraction_table, file.path(csv_name[i], paste0("fraction_", csv_name[i], "_outliers_per_feature.csv")))
    
    # outlier enrichment analysis
    outlier_analysis_out <- outlier_analysis(grouptablist = count_outliers_out$grouptablist,
                                             fraction_table = fraction_table,
                                             fraction_samples_cutoff = fraction_samples_cutoff,
                                             write_out_tables = TRUE,
                                             outfilepath = file.path(csv_name[i], "/"))
    
    results[[i]] = list(fraction_table = fraction_table, outlier_analysis_out = outlier_analysis_out)
    print(paste0(csv_name[i], " outlier analysis complete"))
  }
  return(results)
}

generate_heatmap_annotations = function(gct, heatmap_annotations_columns, annotations_columns_of_interest){

  # create and order annotations for heatmap display
  annotations_columns = strsplit(heatmap_annotations_columns, ",")[[1]]
  
  annotations_ordered = gct@cdesc %>%
    select(all_of(annotations_columns)) %>%
    rownames_to_column("rowname") %>%
    arrange_at(annotations_columns_of_interest) %>%
    column_to_rownames("rowname")
  
  return(annotations_ordered)
}

generate_outliers_heatmaps = function(annotations, outliers_results_pos_neg, fdrcutoffvalue, heatmap_annotations){
  
  csv_name = c("negative", "positive")
  heatmap_names = names(data.frame(annotations))
  
  for (pos_neg in 1:length(outliers_results_pos_neg)){
    deva_output = outliers_results_pos_neg[[pos_neg]]
    outlier_analysis_out = deva_output$outlier_analysis_out
    heatmaplist = NULL
    
    # adapted from blacksheepr outlier_heatmap function
    for (i in 1:length(outlier_analysis_out)){
      intable = data.frame(outlier_analysis_out[[i]])
      intable[intable == ""] = NA
      
      # select column that contains fdr values for the "in" group, get rownames (genes) that meet fdr cutoff
      fdrcols = grep("fdr", colnames(intable), value = TRUE)
      fdrcols = fdrcols[which(str_detect(fdrcols, "not", negate = TRUE))]
      
      # genes of interest for heatmap
      GOI = as.character(intable[which(intable[,fdrcols]<fdrcutoffvalue),1])
      title = paste0("Significant genes in ", csv_name[pos_neg], " outlier analysis for ", heatmap_names[i], ", FDR cutoff value = ", fdrcutoffvalue)
      GOI_list = append(title, GOI)
      
      # save list of genes of interest
      write.table(GOI_list, file.path(csv_name[pos_neg], paste0("significant_genes_ ", csv_name[pos_neg], "_outlier_analysis_", heatmap_names[i], ".txt")), quote = FALSE, row.names = FALSE, col.names = FALSE)

      if(length(GOI) > 0) {
        subset_fraction_table = deva_output$fraction_table[GOI,
                                                            rownames(heatmap_annotations), drop=FALSE]
        
        col_heatmap_annotations = annotationlist_builder(heatmap_annotations)
        
        heatmap = create_heatmap(counttab = subset_fraction_table,
                                 colmetatable = heatmap_annotations,
                                 colannotationlist = col_heatmap_annotations,
                                 colclusterparam = FALSE,
                                 rowclusterparam = FALSE,
                                 nameparam = paste0(csv_name[pos_neg], " outlier analysis for ", heatmap_names[i], ", FDR cutoff value = ", fdrcutoffvalue))
        
        pdf(file.path(csv_name[pos_neg], paste0(csv_name[pos_neg], "_outlier_analysis_", heatmap_names[i], ".pdf")))
        print(heatmap)
        dev.off()
        
        print(paste0(csv_name[pos_neg], " outlier analysis heatmap for ", heatmap_names[i], " complete"))
        
      }
    }
  }
}

run_deva_all = function(annotations, data_values, fraction_samples_cutoff,
                        gct, groups_file, fdrcutoffvalue){
  
  # make annotation groupings
  groupings = comparison_groupings(annotations)
  
  # count outliers and run enrichment analysis for both positive and negative outliers
  outliers_results_pos_neg = outlier_count_and_analysis(groupings, data_values, fraction_samples_cutoff)
  
  # generate heatmap annotations
  heatmap_annotations = generate_heatmap_annotations(gct, heatmap_annotations_columns, annotations_columns_of_interest)
  
  # generate heatmaps
  generate_outliers_heatmaps(annotations, outliers_results_pos_neg, fdrcutoffvalue, heatmap_annotations)
  
  return(outliers_results_pos_neg)
    
}
  
run_blacksheep_analysis = function(gct_path, 
                                   GeneSymbol_column, 
                                   identifiers_file, 
                                   groups_file,
                                   fraction_samples_cutoff, 
                                   fdrcutoffvalue){
  
  dir.create("blacksheep")
  setwd("blacksheep")
  
  # read gct file
  gct = parse.gctx(gct_path)

  if (groups_file == NULL) {

    data_values = data.frame(gct@mat)

    # generate and save outlier tables for positive and negative outliers
    blacksheep_out = outlier_count_only(data_values)

  } else {

    #format data and annotations
    data_values = create_values_input(gct, GeneSymbol_column, identifiers_file)
    annotations = create_annotations_input(gct, groups_file)

    # run deva for positive and negative outliers, save tables and generate heatmaps
    blacksheep_out = run_deva_all(annotations, data_values, fraction_samples_cutoff,
                          gct, groups_file, fdrcutoffvalue)
  
  }

  return(blacksheep_out)
}

blksheep = run_blacksheep_analysis(gct_path = gct_path,
                                   GeneSymbol_column = GeneSymbol_column,
                                   identifiers_file = identifiers_file,
                                   groups_file = groups_file,
                                   fraction_samples_cutoff = fraction_samples_cutoff,
                                   fdrcutoffvalue = fdrcutoffvalue)
