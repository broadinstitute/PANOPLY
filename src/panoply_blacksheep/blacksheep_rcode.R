#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
library(blacksheepr)
library(cmapR)
library(dplyr)
library(tibble)
library(stringr)
library(yaml)
library(RColorBrewer)
library(grid)

args <- commandArgs(TRUE)

gct_path = as.character(args[1])
yaml_file = as.character(args[2])

yaml_params = read_yaml(yaml_file)
GeneSymbol_column = yaml_params$global_parameters$gene_mapping$gene_id_col
SampleID_column = yaml_params$DEV_sample_annotation$sample_id_col_name
apply_identifiers_filter = yaml_params$panoply_blacksheep$apply_filtering
identifiers_file = yaml_params$panoply_blacksheep$identifiers_file
groups_file = yaml_params$panoply_blacksheep$groups_file
fraction_samples_cutoff = yaml_params$panoply_blacksheep$fraction_samples_cutoff
fdrcutoffvalue = yaml_params$panoply_blacksheep$fdr_value

create_values_input = function(gct, GeneSymbol_column, apply_identifiers_filter, identifiers_file){
  
  # extract data table and replace protein name with gene symbol for downstream aggregation and analysis
  data_values = data.frame(gct@mat) %>%
    rownames_to_column("rowname")

  genesymbol = gct@rdesc %>%
    rownames_to_column("rowname") %>%
    select(all_of(GeneSymbol_column), rowname) %>%
    left_join(data_values, by = "rowname")

  if (isTRUE(apply_identifiers_filter)){
    if (is.null(identifiers_file)){
      identifiers_file = "/prot/proteomics/Projects/PGDAC/src/kinase_list.txt"
      print("Kinases filter applied")
    } else {
      print("Applying user-supplied file filter")
    }
    identifiers = read.delim(identifiers_file, header = FALSE, sep="\t")
    identifiers = as.character(identifiers$V1)
    genesymbol = genesymbol %>%
      filter(genesymbol[,GeneSymbol_column] %in% identifiers)
  } else {
    print("No filter applied")
  }
  
  genesymbol[,GeneSymbol_column] = paste(genesymbol[,GeneSymbol_column], row.names(genesymbol), sep = "-")
  
  genesymbol_values = genesymbol %>%
    column_to_rownames(GeneSymbol_column) %>%
    select(-rowname)
  
  return(genesymbol_values)
}

create_groupsfile_annotations = function(gct, groups_file, SampleID_column, data_values){
  
  annotations = read.csv(groups_file, stringsAsFactors = F, na.strings=c("", "NA")) %>%
    column_to_rownames(SampleID_column)
  annotations[is.na(annotations)] = "NA"
  
  annotations = annotations[names(data_values),]
  
  return(annotations)
}

create_binary_annotations = function(annotations){
  
  binary_annotations = make_comparison_columns(annotations)
  
  cols = vector()
  for (i in colnames(binary_annotations)){
    if (1 %in% table(binary_annotations[,i])){
      cols = append(i, cols)
    }
  }
  
  binary_annotations_revised = binary_annotations[, !(colnames(binary_annotations) %in% cols)]
  
  return(binary_annotations_revised)
  
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
  }

  return(outlier_table)
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
    write.csv(fraction_table, file.path(csv_name[i], paste0("fraction_", csv_name[i], "_outliers_per_genesymbol.csv")))
    
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

check_annot_colors = function(color_list, annotations, blank.color = "#FFFFFF"){
  # color_list: colors from yaml parsed as list
  # annotations: groups file annotations
  # blank.color: color for blank values, default = white
  ## note: blank values in groups annotations should already be converted to "NA" by create_groupsfile_annotations function
  
  keep.idx <- names(color_list) %in% colnames(annotations)
  cdesc.color <- color_list[ keep.idx ]
  
  for(i in 1:length(cdesc.color)){
    
    ## cdesc column
    cdesc.column <- names(cdesc.color)[i]
    
    cdesc.levels <- annotations[, cdesc.column] %>% unique
    
    #replace any blank or NA values to character "NA"
    cdesc.levels[ nchar(as.character(cdesc.levels)) == 0 | is.na(cdesc.levels) ] <- "NA"
    
    ## extract colors specified for levels of cdesc.column
    color.tmp <- color_list[[i]]
    
    ## check whether all levels were assigned
    ## fill with random colors
    if(sum(!cdesc.levels %in% names(color.tmp)) > 0){
      missed_names = cdesc.levels[which(!cdesc.levels %in% names(color.tmp))]
      
      if("NA" %in% missed_names){
        color.tmp["NA"] <- blank.color
        missed_names = missed_names[-which(missed_names == "NA")]
      }
      
      if(length(missed_names)>0){
        missed_colors = brewer.pal(length(missed_names), "Set3")[1:length(missed_names)]
        names(missed_colors) = missed_names
        
        color.tmp = append(color.tmp, missed_colors)
      }

      cdesc.color[[i]] <- color.tmp
    }
  }
  
  return(cdesc.color)
}

generate_outliers_heatmaps = function(binary_annotations, outliers_results_pos_neg, fdrcutoffvalue, annotations){
  
  csv_name = c("negative", "positive")
  category_names = names(data.frame(binary_annotations))
  
  # parse colors from yaml file for heatmap annotation label
  color = lapply(yaml_params$groups.colors, unlist)
  color2 = check_annot_colors(color_list = color, annotations = annotations)
  
  for (pos_neg in 1:length(outliers_results_pos_neg)){
    deva_output = outliers_results_pos_neg[[pos_neg]]
    outlier_analysis_out = deva_output$outlier_analysis_out
    heatmaplist = NULL
    
    # adapted from blacksheepr outlier_heatmap function
    for (i in 1:length(outlier_analysis_out)){
      
      for (name in category_names){
        if(grepl(name, names(outlier_analysis_out[i]))){
          heatmap_name = name
        }
      }
      
      intable = data.frame(outlier_analysis_out[[i]])
      intable[intable == ""] = NA
      
      # select column that contains fdr values for the "in" group, get rownames (genes) that meet fdr cutoff
      fdrcols = grep("fdr", colnames(intable), value = TRUE)
      fdrcols = fdrcols[which(str_detect(fdrcols, "not", negate = TRUE))]
      intable[,fdrcols] = as.numeric(intable[,fdrcols])
      
      # genes of interest for heatmap
      GOI = as.character(intable[which(intable[,fdrcols]<fdrcutoffvalue),1])
      title = paste0("Significant genes in ", csv_name[pos_neg], " outlier analysis for ", heatmap_name, ", FDR cutoff value = ", fdrcutoffvalue)
      GOI_list = append(title, GOI)
      
      # save list of genes of interest
      write.table(GOI_list, file.path(csv_name[pos_neg], paste0("significant_genes_", csv_name[pos_neg], "_outlier_analysis_", heatmap_name, ".txt")), quote = FALSE, row.names = FALSE, col.names = FALSE)

      if(length(GOI) > 0) {
        
        #order annotations for heatmap
        column = gsub("_[^_]+$", "", heatmap_name)
        annotations_ordered = annotations %>%
          rownames_to_column("rowname") %>%
          arrange_at(column) %>%
          column_to_rownames("rowname")
        
        subset_fraction_table = deva_output$fraction_table[GOI,
                                                            rownames(annotations_ordered), drop=FALSE]
        
        col_heatmap_annotations = annotationlist_builder(annotations_ordered, customcolorlist = color2)
        
        heatmap = create_heatmap(counttab = subset_fraction_table,
                                 colmetatable = annotations_ordered,
                                 colannotationlist = col_heatmap_annotations,
                                 colclusterparam = FALSE,
                                 rowclusterparam = FALSE,
                                 nameparam = paste0(csv_name[pos_neg], " outlier analysis for ", heatmap_name, ", FDR cutoff value = ", fdrcutoffvalue))
        if (dim(annotations_ordered)[1] > 11){
          pdf(file.path(csv_name[pos_neg], paste0(csv_name[pos_neg], "_outlier_analysis_", heatmap_name, ".pdf")), 
              height = unit(10, "cm"),
              width = unit(10, "cm"))
          print(heatmap)
          dev.off()
          
          png(file.path(csv_name[pos_neg], paste0(csv_name[pos_neg], "_outlier_analysis_", heatmap_name, ".png")), 
              width = 3000,
              height = 3000,
              res = 300)
          print(heatmap)
          dev.off()
        } else {
          pdf(file.path(csv_name[pos_neg], paste0(csv_name[pos_neg], "_outlier_analysis_", heatmap_name, ".pdf")))
          print(heatmap)
          dev.off()
          
          png(file.path(csv_name[pos_neg], paste0(csv_name[pos_neg], "_outlier_analysis_", heatmap_name, ".png")), width = 2200, height = 2200, res = 300)
          print(heatmap)
          dev.off()
        }
        
        print(paste0(csv_name[pos_neg], " outlier analysis heatmap for ", heatmap_name, " complete"))
        
      }
    }
  }
}

run_deva_all = function(binary_annotations, data_values, fraction_samples_cutoff, fdrcutoffvalue, annotations){
  
  # make annotation groupings
  groupings = comparison_groupings(binary_annotations)
  
  # count outliers and run enrichment analysis for both positive and negative outliers
  outliers_results_pos_neg = outlier_count_and_analysis(groupings, data_values, fraction_samples_cutoff)
    
  # generate heatmaps
  generate_outliers_heatmaps(binary_annotations, outliers_results_pos_neg, fdrcutoffvalue, annotations)
  
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

  if (is.null(groups_file)) {

    data_values = data.frame(gct@mat)

    # generate and save outlier tables for positive and negative outliers
    blacksheep_out = outlier_count_only(data_values)

  } else {

    #format data and annotations
    data_values = create_values_input(gct, GeneSymbol_column, apply_identifiers_filter, identifiers_file)
    heatmap_annotations = create_groupsfile_annotations(gct, groups_file, SampleID_column, data_values)

    # binarize annotations
    binary_annotations = create_binary_annotations(heatmap_annotations)

    # run deva for positive and negative outliers, save tables and generate heatmaps
    blacksheep_out = run_deva_all(binary_annotations, data_values, fraction_samples_cutoff,
                          fdrcutoffvalue, heatmap_annotations)
  
  }

  return(blacksheep_out)
}

blksheep = run_blacksheep_analysis(gct_path = gct_path,
                                   GeneSymbol_column = GeneSymbol_column,
                                   identifiers_file = identifiers_file,
                                   groups_file = groups_file,
                                   fraction_samples_cutoff = fraction_samples_cutoff,
                                   fdrcutoffvalue = fdrcutoffvalue)
