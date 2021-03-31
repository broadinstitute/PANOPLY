library(rmimp)
library(seqinr)
library(dplyr)
library(tidyr)
library(tibble)
library(cmapR)
library(readr)
library(S4Vectors)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(stringr)

#args <- commandArgs(TRUE)

# mut_data_path = as.character(args[1])
# phospho_path = as.character(args[2])
# fasta_path = as.character(args[3])
# ids_path = as.character(args[4])
# yaml_file = as.character(args[5])

# yaml_params = read_yaml(yaml_file)
# ids = whatever the name of the ids file within the .Rdata file is. default ids?

# phospho file column names
SpectrumMill = TRUE
phosphosite_col = "variableSites"
accession_number_col = "accession_number"

## mutation file column names - should this be required or in yaml?
protein_change_colname = "Protein_Change"
mutation_type_col = "Variant_Classification"
patient_id_col = "Tumor_Sample_Barcode"
RefSeq_col = "refseq_mrna_id"

setwd("C:/Users/karen/mimp")
mut_data_path = "S:/CPTAC3/PGDAC/brca/prospective/v5.4-public/data-freeze/prosp-brca-v5.4-public-BRCA-freeze-v5.final_analysis_set.maf.txt"
phospho_path = "S:/CPTAC3/PGDAC/brca/prospective/v5.4-public/data-freeze/prosp-brca-v5.4-public-phosphoproteome-ratio-norm-NArm.gct"
ids_path = "P:/LabMembers/Abhijeet Mavi/MIMP_final/ids.RData"
fasta_path = "S:/CPTAC2/RefSeq_20160914/RefSeq.20160914_Human_ucsc_hg19_customProDBnr_mito_150contams.fasta"

## phospho file column names - if multiple search engines, uncomment and add other options via if
# if (SpectrumMill){
#   phosphosite_col = "variableSites"
#   accession_number_col = "accession_number"
# } 
#   

# this function prepares fasta file for mimp input
format_fasta_file = function(fasta_path){
  # Preprocessing the fasta file in RefSeq format
  fasta = Biostrings::readAAStringSet(fasta_path)
  names(fasta) <- names(fasta) %>% sub('\\..*', '', .) %>% sub('^>','', .)
  keep = grepl("^NP", names(fasta))
  seqdata = as.list( as.character( fasta[keep] ) )
  
  return(seqdata)
}

# this function prepares the mutation maf file before sample-wise filtering
format_mutation_full = function(mut_data_path, ids, mutation_type_col, RefSeq_col){
  mut_maf = read_tsv(mut_data_path)
  mut_data = as.data.frame(mut_maf) 
  mut_data[, mutation_type_col] = as.character(mut_data[, mutation_type_col])
  mut_data2 = mut_data[which(grepl("missense",mut_data[, mutation_type_col], ignore.case = TRUE)),]
  
  # map NP numbers to NM
  mut_data2$NM_revised = sub('\\..*','',mut_data2[, RefSeq_col])
  
  keep.idx = !is.na(mut_data2$NM_revised)
  
  ids2 = ids %>%
    select(tx_name, pro_name)
  names(ids2) = c("NM_revised", "NP_id")
  
  mut_data3 = mut_data2[keep.idx, ] %>%
    left_join(ids2) %>%
    filter(NP_id %in% names(seqdata))
  
  return (mut_data3)
}

# create sample-wise mutation mimp input, executed in run_mimp_samplewise
create_mutation_input = function(mut_sample, loop_id, protein_change_colname, sample_input_path){
  # Creating the dataframe for mutation file in MIMP format from respective columns of mutation file
  df_mut = mut_sample %>%
    select(NP_id,
           all_of(protein_change_colname))
  names(df_mut)[2] = "mutation"
  df_mut = df_mut %>%  
    mutate(mutation = gsub('p\\.','', mutation))
  
  #Validating if df_mut is in correct format. Eg: NP_027626728
  mut_rows <- which(grepl('^[A-Z][0-9]+[A-Z]$',df_mut$mutation) &
                      (!is.na(df_mut$NP_id)) & (df_mut$NP_id!=''))
  
  # This will be used as an input into MIMP algorithm. Eg: NP_0976666 K226R
  df_mut2 <- df_mut[ mut_rows,]
  write.csv(df_mut, file.path(sample_input_path, paste0(i, "_mutation_mimp_input.csv")), row.names = FALSE)
  
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

# create sample-wise phospho mimp input, executed in run_mimp_samplewise
create_phospho_input = function(phospho_mat, phospho_rdesc, loop_id, phosphosite_col, path){
  # Creating the phospho data in MIMP format
  phos_mat_id = !is.na(phospho_mat[,loop_id])
  
  df_phospho = phospho_rdesc[phos_mat_id,] %>%
    select(NP_revised,
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

run_mimp_samplewise = function(phospho_cid, mut_data, seqdata, phos_rdesc, phos_mat2, 
                               phosphosite_col, protein_change_colname){

  # set up empty lists/data frames to save or concatenate loop output
  log_no_mut = "Samples with no mutations (MIMP not run): \n"
  log_no_mimp_res = "Samples with no predicted kinase rewiring results from MIMP: \n"
  
  full_results = data.frame(gene = as.character(),
                            wt = as.character(),
                            mt = as.character(),
                            mut = as.character(),
                            psite_pos	= as.numeric(),
                            mut_dist = as.numeric(),
                            score_wt = as.numeric(),
                            score_mt = as.numeric(),
                            log_ratio	= as.numeric(),
                            pwm = as.character(),
                            pwm_fam = as.character(),
                            nseqs = as.numeric(),
                            prob = as.numeric(),
                            effect = as.character(),
                            patient_id = as.character())
  
  gain_STY_all = data.frame(NP_id = as.character(),
                             mutation = as.character(),
                             patient_id = as.character())
  lose_STY_all = data.frame(NP_id = as.character(),
                             mutation = as.character(),
                             patient_id = as.character())
  
  # Algorithm to run MIMP on whole dataset, patient-wise
  for (i in phosho_cid){
    
    dir.create(file.path("results_by_sample", i))
    
    sample_input_path = file.path("results_by_sample", i, "mimp_input")
    sample_results_path = file.path("results_by_sample", i, "mimp_results")
    sample_mutinfo_path = file.path("results_by_sample", i, "mutation_info")
    dir.create(sample_input_path)
    dir.create(sample_results_path)
    dir.create(sample_mutinfo_path)
    
    # select patient-specific mutation entries
    mut_i = mut_data[which(grepl(i,mut_data[, patient_id_col])),]
    write.csv(mut_i, file.path(sample_mutinfo_path, paste0(i, "_mutation_data_all.csv")), row.names = FALSE)
    
    if(dim(mut_i)[1] > 0){
      
      df_mut = create_mutation_input(mut_i, i, protein_change_colname, sample_input_path)
      
      gain_STY = central_STY_change(df_mut, i, "gain", sample_results_path)     
      gain_STY_all = rbind(gain_STY_all, gain_STY)
      
      lose_STY = central_STY_change(df_mut, i, "loss", sample_results_path)
      lose_STY_all = rbind(lose_STY_all, lose_STY)
      
      # save list of mutation reference AAs that don't match with fasta sequence
      AA_mismatch = data.frame(protein = as.character(),
                               mutation = as.character(),
                               mutation_position = as.numeric(),
                               mutation_reference_AA = as.character(),
                               fasta_AA = as.character())
      # AA_windows = data.frame(protein = as.character(),
      #                         mutation = as.character(),
      #                         window_containing_STY = as.character())
      for (j in 1:length(df_mut$NP_id)){
        AA_pos = as.numeric(gsub("[A-Z]", "", df_mut$mutation[j]))
        if (str_sub(seqdata[df_mut$NP_id[j]], AA_pos, AA_pos) != str_sub(df_mut$mutation[j], 1, 1)){
          j_df = data.frame(protein = df_mut$NP_id[j],
                            mutation = df_mut$mutation[j],
                            AA_position = AA_pos,
                            AA_mutation_reference = str_sub(df_mut$mutation[j], 1, 1),
                            AA_fasta = str_sub(seqdata[df_mut$NP_id[j]], AA_pos, AA_pos))
          AA_mismatch = rbind(AA_mismatch, j_df)
        }
        # AA_window = str_sub(seqdata[df_mut2$NP_id[j]], AA_pos - 7, AA_pos + 7)
        # if (grepl("S|T|Y", AA_window)){
        #   new_df = data.frame(protein = df_mut2$NP_id[j],
        #                       mutation = df_mut2$mutation[j],
        #                       window_containing_STY = AA_window)
        #   AA_windows = rbind(AA_windows, new_df)
        # }
      }
      
      if (dim(AA_mismatch)[1] > 0){
        write.csv(AA_mismatch, file.path(sample_mutinfo_path, paste0(i, "_mismatchedAA_mutation_vs_fasta.csv")), row.names = FALSE)
      }
      
      # if (dim(AA_windows)[1] > 0){
      #   write.csv(AA_windows, file.path(sample_mutinfo_path, paste0(i, "_mutations_with_neighboring_STY.csv")), row.names = FALSE)
      # }
      
      df_phospho = create_phospho_input(phos_mat2, phos_rdesc, i, phosphosite_col, sample_input_path)

      # note any protein identifiers from phospho data that are not in the fasta file
      prot_not_in_fasta = list(df_phospho$NP_revised[which(!(df_phospho$NP_revised %in% names(seqdata)))])
      if (length(prot_not_in_fasta[[1]]) > 0){
        names(prot_not_in_fasta) = "phospho identifiers without fasta match"
        write.csv(prot_not_in_fasta, file.path(sample_mutinfo_path, paste0(i, "_phospho_identifiers_not_in_fasta.csv")), row.names = FALSE)
      }
      
      # Run MIMP (patient-specific)
      results_i = mimp(muts = df_mut, seqs = seqdata, psites = df_phospho)
      
      if (!is.null(results_i)){
        
        # generate HTML report
        if (file.exists(file.path(BASE_DIR, 'html', 'MIMP_results.html'))){
          file.copy(file.path(BASE_DIR, 'html', 'MIMP_results.html'), paste0(i, "_mimp_output_kinase_rewiring_events.html"))
        }
        
        results_i = results_i %>%
          distinct() %>%
          mutate(patient_id = i)
        write.csv(results_i, file = file.path(sample_results_path, paste(i, "mimp_output_kinase_rewiring_events.csv", sep = "_")), row.names = FALSE)
        full_results = rbind(full_results, results_i)
        
      } else {
        log_no_mimp_res = paste(log_no_mimp_res, i, sep = "\n")
      }
    } else {
      log_no_mut = paste(log_no_mut, i, sep = "\n")
    }
    
  }
  
  if(dim(full_results)[1]>0){
    write.csv(full_results, "mimp_output_kinase_rewiring_events_all.csv", row.names = FALSE)
  } else{
    print("No predicted kinase rewiring events in any samples")
  }
  
  if(dim(gain_site_all)[1]>0){
    write.csv(gain_site_all, "all_central_residue_STY_gain.csv", row.names = FALSE)
  } else{
    print("No mutations in any sample that cause gain of new phosphosite")
  }
  
  if(dim(loss_site_all)[1]>0){
    write.csv(loss_site_all, "all_central_residue_STY_loss.csv", row.names = FALSE)
  } else{
    print("No mutations in any sample that cause loss of a phosphosite")
  }
  
  writeLines(log_no_mut, con = "samples_without_mutations.txt")
  writeLines(log_no_mimp_res, con = "samples_without_kinase_rewiring_results.txt")
  
  return(full_results)
}  

mimp_heatmap_function = function(df, kinase_column){

  heatmap_df = df %>%
    select(all_of(kinase_column),log_ratio, patient_id) %>%
    group_by_at(c(kinase_column, "patient_id")) %>%
    slice(which.max(abs(log_ratio))) %>%
    spread_(key = kinase_column, value = "log_ratio") %>%
    column_to_rownames("patient_id") %>%
    t()
  write.csv(heatmap_df, paste("kinase_rewiring_events_matrix", kinase_column, "level.csv"))

  results.max = max(df$log_ratio, na.rm = TRUE)
  results.min = min(df$log_ratio, na.rm = TRUE)
  
  color = colorRamp2(c(results.min, 0, results.max), c("blue", "white", "red"))
  
  heatmap_all = Heatmap(heatmap_df, cluster_rows = FALSE, cluster_columns = FALSE, 
                        col = color, column_title = "Sample", row_title = "Kinase",
                        column_title_side = "bottom", column_title_gp = gpar(fontsize = 10),
                        row_title_side = "left", row_title_gp = gpar(fontsize = 10),
                        row_names_side = "left", row_names_gp = gpar(fontsize=6),
                        column_names_gp = gpar(fontsize = 6), 
                        heatmap_legend_param = list(title = "log_ratio"),
                        width = ncol(full_results_edit)*unit(2, "mm"),
                        height = nrow(full_results_edit)*unit(2, "mm"))
  
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

run_mimp = function(fasta_path, phospho_path, SpectrumMill, accession_number_col,
                    mut_data_path, ids, mutation_type_col, RefSeq_col, 
                    phosphosite_col, protein_change_colname){

  dir.create("results_dir")
  setwd("results_dir")
  dir.create("results_by_sample")

  # Loading NP data, different for different datatypes
  load(ids_path)
  
  # prepare fasta file for mimp input
  seqdata = format_fasta_file(fasta_path)

  # read phospho gct file and remove Xs in front of patient IDs
  phospho_gct = parse.gctx(phospho_path)
  phos_rdesc <- phospho_gct@rdesc
  phos_mat <- phospho_gct@mat
  colnames(phos_mat) <- sub('^X','', colnames(phos_mat))
  phos_cid <- phospho_gct@cid %>% sub('^X', '', .)
  
  # filter for fully localized sites if SpectrumMill (# localized = # actual sites)
  if (SpectrumMill){
    phos_rdesc = phos_rdesc %>%
      filter(Best_numActualVMSites_sty == Best_numLocalizedVMsites_sty)
    rownames(phos_rdesc) = phos_rdesc$id
  }
  
  phos_mat2 = phos_mat[rownames(phos_rdesc),]
  phos_rdesc$NP_revised = sub('\\..*','',phos_rdesc[, accession_number_col])
  
  # prepare mutation maf file before for loop sample-wise processing
  mut_data = format_mutation_full(mut_data_path, ids, mutation_type_col, RefSeq_col)

  # prepare sample-wise mutation and phospho dfs, and then run mimp on each sample
  heatmap_df = run_mimp_samplewise(phos_cid, mut_data, seqdata, phos_rdesc, phos_mat2,
                                   phosphosite_col, protein_change_colname)

  # create heatmap of predicted kinases results
  generate_mimp_heatmap(heatmap_df)
  
  return(heatmap_df)
}

mimp_results = run_mimp(fasta_path = fasta_path, 
                        phospho_path = phospho_path,
                        SpectrumMill = SpectrumMill, 
                        accession_number_col = accession_number_col,
                        mut_data_path = mut_data_path, 
                        ids = ids,
                        mutation_type_col = mutation_type_col, 
                        RefSeq_col = RefSeq_col,
                        phosphosite_col = phosphosite_col, 
                        protein_change_colname = protein_change_colname)