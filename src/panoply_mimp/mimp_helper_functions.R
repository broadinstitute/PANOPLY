## functions used in panoply_mimp.R

# this function prepares fasta file for mimp input
format_fasta_file = function(fasta_path){
  # Preprocessing the fasta file in RefSeq format
  fasta = Biostrings::readAAStringSet(fasta_path)
  names(fasta) <- names(fasta) %>% sub('\\..*', '', .) %>% sub('^>','', .)
  keep = grepl("^NP", names(fasta))
  seqdata = as.list( as.character( fasta[keep] ) )
  
  return(seqdata)
}

# this function formats the phospho GCT row metadata
format_phospho_rdesc = function(gct, search_engine, accession_number_col){
  phos_rdesc = gct@rdesc
  
  # filter for fully localized sites if SpectrumMill (# localized = # actual sites)
  if (search_engine == "SpectrumMill"){
    phos_rdesc = phos_rdesc %>%
      filter(Best_numActualVMSites_sty == Best_numLocalizedVMsites_sty)
    rownames(phos_rdesc) = phos_rdesc$id
  }
  
  phos_rdesc$protein_id = sub('\\..*','',phos_rdesc[, accession_number_col])
  
  return(phos_rdesc)
}

# this function formats the phopsho data matrix
format_phospho_mat = function(gct, rdesc){
  phos_mat = gct@mat
  colnames(phos_mat) <- sub('^X','', colnames(phos_mat))
  phos_mat = phos_mat[rownames(rdesc),]
  
  return(phos_mat)
}

# this function prepares the mutation maf file before sample-wise filtering
format_mutation_full = function(mut_data_path, ids, mutation_type_col, RefSeq_col, seqdata){
  mut_maf = read_tsv(mut_data_path)
  mut_data = as.data.frame(mut_maf) 
  mut_data[, mutation_type_col] = as.character(mut_data[, mutation_type_col])
  mut_data2 = mut_data[which(grepl("missense",mut_data[, mutation_type_col], ignore.case = TRUE)),]
  
  # map NP numbers to NM
  mut_data2$transcript_id = sub('\\..*','',mut_data2[, RefSeq_col])
  
  keep.idx = !is.na(mut_data2$transcript_id)
  
  ids2 = ids %>%
    select(tx_name, pro_name)
  names(ids2) = c("transcript_id", "protein_id")
  
  mut_data3 = mut_data2[keep.idx, ] %>%
    left_join(ids2) %>%
    filter(protein_id %in% names(seqdata))
  
  return (mut_data3)
}

# create sample-wise mutation mimp input, executed in run_mimp_samplewise
create_mutation_input = function(mut_sample, loop_id, protein_change_colname, sample_input_path){
  # Creating the dataframe for mutation file in MIMP format from respective columns of mutation file
  df_mut = mut_sample %>%
    select(protein_id,
           all_of(protein_change_colname))
  names(df_mut)[2] = "mutation"
  df_mut = df_mut %>%  
    mutate(mutation = gsub('p\\.','', mutation))
  
  #Validating if df_mut is in correct format. Eg: NP_027626728
  mut_rows <- which(grepl('^[A-Z][0-9]+[A-Z]$',df_mut$mutation) &
                      (!is.na(df_mut$protein_id)) & (df_mut$protein_id!=''))
  
  # This will be used as an input into MIMP algorithm. Eg: NP_0976666 K226R
  df_mut2 <- df_mut[ mut_rows,]
  write.csv(df_mut, file.path(sample_input_path, paste0(loop_id, "_mutation_mimp_input.csv")), row.names = FALSE)
  
  return(df_mut2)
}

# find mutations that directly remove or add an STY, executed in run_mimp_samplewise
central_STY_change = function(mimp_mut_input, loop_id, gain_or_loss, sample_results_path){
  if (gain_or_loss == "gain"){
    STY_change = mimp_mut_input %>%
      filter(!(str_sub(mutation, 1, 1) %in% c("S", "T", "Y"))) %>%
      filter(str_sub(mutation, -1) %in% c("S", "T", "Y"))
  } else if (gain_or_loss == "loss"){
    STY_change = mimp_mut_input %>%
      filter(!(str_sub(mutation, -1) %in% c("S", "T", "Y"))) %>%
      filter(str_sub(mutation, 1, 1) %in% c("S", "T", "Y"))
  }
  
  if (dim(STY_change)[1] > 0){
    write.csv(STY_change, file.path(sample_results_path, paste0(loop_id, "_central_residue_STY_", gain_or_loss, ".csv")), row.names = FALSE)
  }
  
  STY_change = STY_change %>% mutate(patient_id = loop_id)
  
  return(STY_change)
}

# function to save list of mutation reference AAs that don't match with fasta sequence
save_mismatch_AA = function(df_mut, seqdata, sample_mutinfo_path, loop_id){
  AA_mismatch = data.frame(protein_id = as.character(),
                           mutation = as.character(),
                           mutation_AA_position = as.numeric(),
                           mutation_reference_AA = as.character(),
                           fasta_AA = as.character())
  
  for (j in 1:length(df_mut$protein_id)){
    AA_pos = as.numeric(gsub("[A-Z]", "", df_mut$mutation[j]))
    if (str_sub(seqdata[df_mut$protein_id[j]], AA_pos, AA_pos) != str_sub(df_mut$mutation[j], 1, 1)){
      j_df = data.frame(protein_id = df_mut$protein_id[j],
                        mutation = df_mut$mutation[j],
                        mutation_AA_position = AA_pos,
                        mutation_reference_AA = str_sub(df_mut$mutation[j], 1, 1),
                        fasta_AA = str_sub(seqdata[df_mut$protein_id[j]], AA_pos, AA_pos))
      AA_mismatch = rbind(AA_mismatch, j_df)
    }
  }
  
  if (dim(AA_mismatch)[1] > 0){
    write.csv(AA_mismatch, file.path(sample_mutinfo_path, paste0(loop_id, "_mismatchedAA_mutation_vs_fasta.csv")), row.names = FALSE)
  }
  
  return(AA_mismatch)
}

# create sample-wise phospho mimp input, executed in run_mimp_samplewise
create_phospho_input = function(phospho_mat, phospho_rdesc, loop_id, phosphosite_col, path){
  # Creating the phospho data in MIMP format
  phos_mat_id = !is.na(phospho_mat[,loop_id])
  
  df_phospho = phospho_rdesc[phos_mat_id,] %>%
    select(protein_id,
           all_of(phosphosite_col))
  names(df_phospho)[2] = "phosphosite"
  df_phospho = df_phospho %>%
    mutate(phosphosite = gsub(' $','', phosphosite)) %>%
    mutate(phosphosite = strsplit(phosphosite, " ")) %>%
    unnest(phosphosite) %>%
    filter(grepl("S|T|Y|s|t|y", phosphosite)) %>%
    mutate(phosphosite = gsub("S|T|Y|s|t|y", "", phosphosite))
  write.csv(df_phospho, file.path(path, paste0(loop_id, "_phospho_mimp_input.csv")), row.names = FALSE)
  
  return(df_phospho)
}

# function to save list of mutations that fall in +/- 7 AA window of phosphosite (pSNV)
muts_in_phosphosite_window = function(df_phospho, df_mut, loop_id, AA_mismatch, sample_mutinfo_path){
  pSNV = data.frame(protein_id = as.character(),
                    mutation = as.character(),
                    phosphosite_position = as.numeric(),
                    mut_distance_from_phosphosite = as.numeric())
  
  for (j in 1:dim(df_mut)[1]){
    mut_position = as.numeric(gsub("[A-Z]", "", df_mut$mutation[j]))
    sites = df_phospho[which(df_phospho$protein_id == df_mut$protein_id[j]),]
    if (dim(sites)[1]>0){
      for (k in sites$phosphosite){
        distance = abs(as.numeric(k) - mut_position)
        if (distance <= 7){
          add_pSNV = data.frame(protein_id = df_mut$protein_id[j],
                                mutation = df_mut$mutation[j],
                                phosphosite_position = as.numeric(k),
                                mut_distance_from_phosphosite = distance)
          pSNV = rbind(pSNV, add_pSNV)
        }
      }
    }
  }
  
  if (dim(pSNV)[1] > 0){
    pSNV = pSNV %>% mutate(patient_id = loop_id)
    write.csv(pSNV, file.path(sample_mutinfo_path, paste0(loop_id, "_pSNV_mutations_within_7AA_of_phosphosite.csv")), row.names = FALSE)
    
    if (length(intersect(as.character(pSNV$mutation), as.character(AA_mismatch$mutation))) > 0){
      pSNV2 = pSNV %>% left_join(AA_mismatch)
      write.csv(pSNV2, file.path(sample_mutinfo_path, paste0(loop_id, "_pSNV_mutations_within_7AA_of_phosphosite_including_mismatch_info.csv")), row.names = FALSE)
    }
  }
  
  return(pSNV)
}

mimp_heatmap_function = function(df, kinase_column){
  
  heatmap_df = df %>%
    select(all_of(kinase_column),log_ratio, patient_id) %>%
    group_by_at(c(kinase_column, "patient_id")) %>%
    slice(which.max(abs(log_ratio))) %>%
    spread_(key = kinase_column, value = "log_ratio") %>%
    column_to_rownames("patient_id") %>%
    t()
  write.csv(heatmap_df, paste("kinase_rewiring_events_matrix", kinase_column, "level.csv", sep = "_"))
  
  results.max = max(df$log_ratio, na.rm = TRUE)
  results.min = min(df$log_ratio, na.rm = TRUE)
  
  color = colorRamp2(c(results.min, 0, results.max), c("blue", "white", "red"))
  
  heatmap_all = Heatmap(heatmap_df, cluster_rows = FALSE, cluster_columns = FALSE, 
                        col = color, column_title = "Sample", row_title = "Kinase",
                        column_title_side = "bottom", column_title_gp = gpar(fontsize = 10),
                        row_title_side = "left", row_title_gp = gpar(fontsize = 10),
                        row_names_side = "left", row_names_gp = gpar(fontsize=6),
                        column_names_gp = gpar(fontsize = 6), 
                        heatmap_legend_param = list(title = "log_ratio"))
  
  # print heatmap
  pdf(paste('mimp_results_predicted_kinase_rewiring', kinase_column, 'level_heatmap.pdf', sep = "_"))
  draw(heatmap_all)
  # barplot(t(unmatched_NP),xlab = 'Patient_ID', ylab = 'Percentage of RefSeq IDs matched')
  dev.off()
}

generate_mimp_heatmap = function(full_results){
  
  full_results_edit = full_results %>%
    select(pwm, gene, mut, log_ratio, patient_id) %>%
    rename_(kinase = "pwm") %>%
    filter(!is.na(log_ratio)) %>%
    mutate(kinase_gene_mut = paste(kinase, gene, mut, sep = "_"))
  
  #generate kinase-level heatmap
  mimp_heatmap_function(full_results_edit, "kinase")
  
  #generate kinase + mutation level heatmap
  mimp_heatmap_function(full_results_edit, "kinase_gene_mut")
}