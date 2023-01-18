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
library(yaml)

args <- commandArgs(TRUE)

mut_data_path = as.character(args[1])
phospho_path = as.character(args[2])
fasta_path = as.character(args[3])
ids_path = as.character(args[4])
yaml_file = as.character(args[5])

yaml_params = read_yaml(yaml_file)
groups_file_path = yaml_params$panoply_mimp$groups_file_path
groups_file_SampleID_column = yaml_params$DEV_sample_annotation$sample_id_col_name
search_engine = yaml_params$panoply_mimp$search_engine
phosphosite_col = yaml_params$panoply_mimp$phosphosite_col
protein_id_col = yaml_params$panoply_mimp$protein_id_col
protein_id_type = yaml_params$panoply_mimp$protein_id_type
mutation_AA_change_colname = yaml_params$panoply_mimp$mutation_AA_change_colname
mutation_type_col = yaml_params$panoply_mimp$mutation_type_col
sample_id_col = yaml_params$panoply_mimp$sample_id_col 
transcript_id_col = yaml_params$panoply_mimp$transcript_id_col

## phospho file column names - can add other options as we add new search engines
if (search_engine == "SpectrumMill"){
  phosphosite_col = "variableSites"
  protein_id_col = "accession_number"
  protein_id_type = "ENSEMBLPROT"
} else if (search_engine == "other"){
  if (is.null(phosphosite_col)){
    stop("Please update the 'phosphosite_col' field in the yaml file to the name of the row-metadata field in the phopho GCT that indicates phosphosite information e.g. S18s.")
  }
  if (is.null(protein_id_col)){
    stop("Please update the 'protein_id_col' field in the yaml file to the name of the row-metadata field in the phopho GCT that indicates protein accession number e.g. NP_001611.1.")
  }
  if (is.null(protein_id_type)){
    warning("The 'protein_id_type' field has been left blank, and will be assumed 'REFSEQ' by default. If the 'protein_id_col' column uses a different protein ID type, please update the 'protein_id_type' field in the yaml file.")
    protein_id_type="REFSEQ"
  }
} else {
  stop("Please enter a valid option for 'search_engine' parameter. Options are 'SpectrumMill' or 'other.'")
}

# source helper functions
source("/prot/proteomics/Projects/PGDAC/src/mimp_helper_functions.R")

run_mimp_samplewise = function(phospho_cid, mut_data, seqdata, phos_rdesc, phos_mat, 
                               phosphosite_col, mutation_AA_change_colname){
  
  # set up empty lists/data frames to save or concatenate loop output
  log_no_mut = "Samples with no mutations (MIMP not run): \n"
  log_no_mimp_res = "Samples with no predicted kinase rewiring results from MIMP: \n"
  
  full_results = data.frame(protein_id = as.character(),
                            wild_type_seq = as.character(),
                            mut_seq = as.character(),
                            mutation = as.character(),
                            phosphosite_position = as.numeric(),
                            mut_dist = as.numeric(),
                            score_wt = as.numeric(),
                            score_mt = as.numeric(),
                            log_ratio	= as.numeric(),
                            kinase_pwm = as.character(),
                            kinase_pwm_fam = as.character(),
                            nseqs = as.numeric(),
                            prob = as.numeric(),
                            effect = as.character(),
                            sample_id = as.character())
  
  gain_STY_all = data.frame(protein_id = as.character(),
                            mutation = as.character(),
                            sample_id = as.character())
  lose_STY_all = data.frame(protein_id = as.character(),
                            mutation = as.character(),
                            sample_id = as.character())
  
  AA_mismatch_all = data.frame(protein_id = as.character(),
                               mutation = as.character(),
                               mutation_AA_position = as.numeric(),
                               mutation_reference_AA = as.character(),
                               fasta_AA = as.character())
  
  # Algorithm to run MIMP on whole dataset, sample-wise
  for (i in phospho_cid){
    
    dir.create(file.path("results_by_sample", i))
    
    sample_input_path = file.path("results_by_sample", i, "mimp_input")
    sample_results_path = file.path("results_by_sample", i, "mimp_results")
    sample_mutinfo_path = file.path("results_by_sample", i, "mutation_info")
    dir.create(sample_input_path)
    dir.create(sample_results_path)
    dir.create(sample_mutinfo_path)
    
    # select sample-specific mutation entries
    mut_i = mut_data[which(grepl(i,mut_data[, sample_id_col])),]
    write.csv(mut_i, file.path(sample_mutinfo_path, paste0(i, "_mutation_data_all.csv")), row.names = FALSE)
    
    if(dim(mut_i)[1] > 0){
      
      df_mut = create_mutation_input(mut_i, i, mutation_AA_change_colname, sample_input_path)
      
      gain_STY = STY_change(df_mut, i, "gain", sample_mutinfo_path)     
      gain_STY_all = rbind(gain_STY_all, gain_STY)
      
      lose_STY = STY_change(df_mut, i, "loss", sample_mutinfo_path)
      lose_STY_all = rbind(lose_STY_all, lose_STY)
      
      # save list of mutation reference AAs that don't match with fasta sequence
      AA_mismatch = save_mismatch_AA(df_mut, seqdata, sample_mutinfo_path, i)
      AA_mismatch_all = rbind(AA_mismatch_all, AA_mismatch)
      
      # create phospho input to mimp
      df_phospho = create_phospho_input(phos_mat, phos_rdesc, i, phosphosite_col, sample_input_path)
      
      # save list of mutations that fall in +/- 7 AA window of phosphosite (pSNV)
      pSNV = muts_in_phosphosite_window(df_phospho, df_mut, i, AA_mismatch, sample_mutinfo_path)

      # note any protein identifiers from phospho data that are not in the fasta file
      prot_not_in_fasta = list(df_phospho$protein_id[which(!(df_phospho$protein_id %in% names(seqdata)))])
      if (length(prot_not_in_fasta[[1]]) > 0){
        names(prot_not_in_fasta) = "phospho identifiers without fasta match"
        write.csv(prot_not_in_fasta, file.path(sample_mutinfo_path, paste0(i, "_phospho_identifiers_not_in_fasta.csv")), row.names = FALSE)
      }
      
      # Run MIMP (sample-specific)
      results_i = mimp(muts = df_mut, seqs = seqdata, psites = df_phospho, display.results = TRUE)
      
      if (!is.null(results_i)){
        
        # generate HTML report
        BASE_DIR = system.file("extdata", "", package = "rmimp")
        if (file.exists(file.path(BASE_DIR, 'html', 'MIMP_results.html'))){
          file.copy(file.path(BASE_DIR, 'html'), sample_results_path, recursive = TRUE)
        }
        
        results_i = results_i %>%
          distinct() %>%
          mutate(sample_id = i) %>%
          dplyr::rename(kinase_pwm = pwm,
                        kinase_pwm_fam = pwm_fam,
                        protein_id = gene,
                        wild_type_seq = wt,
                        mut_seq = mt,
                        mutation = mut,
                        phosphosite_position = psite_pos)
        write.csv(results_i, file = file.path(sample_results_path, paste(i, "mimp_output_kinase_rewiring_events.csv", sep = "_")), row.names = FALSE)
        full_results = rbind(full_results, results_i)
        
        # add kinase rewiring events to pSNV CSV file
        for_psnv = results_i %>%
          mutate(kinase_rewiring_events = paste(kinase_pwm, effect, sep = "_")) %>%
          select(protein_id, mutation, phosphosite_position, kinase_rewiring_events) %>%
          group_by(protein_id, mutation, phosphosite_position) %>%
          summarise(num_kinase_rewiring_events = n(),
                    kinase_rewiring_events = paste(kinase_rewiring_events, collapse = ", "))
        
        pSNV = left_join(pSNV, for_psnv) 
        pSNV$num_kinase_rewiring_events[is.na(pSNV$num_kinase_rewiring_events)] = ""
        pSNV$kinase_rewiring_events[is.na(pSNV$kinase_rewiring_events)] = ""
        
        write.csv(pSNV, file.path(sample_mutinfo_path, paste0(i, "_pSNV_mutations_within_7AA_of_phosphosite_including_kinase_rewiring.csv")), row.names = FALSE)
        
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
  
  if(dim(gain_STY_all)[1]>0){
    write.csv(gain_STY_all, "all_STY_gain_mutation.csv", row.names = FALSE)
  } else{
    print("No mutations in any sample that cause gain of new STY")
  }
  
  if(dim(lose_STY_all)[1]>0){
    write.csv(lose_STY_all, "all_STY_loss_mutation.csv", row.names = FALSE)
  } else{
    print("No mutations in any sample that cause loss of STY")
  }
  
  if(dim(AA_mismatch_all)[1]>0){
    write.csv(AA_mismatch_all, "all_mismatchedAA_mutation_vs_fasta.csv", row.names = FALSE)
  } else{
    print("No mismatches between mutation reference and fasta sequence amino acid")
  }
  
  writeLines(log_no_mut, con = "samples_without_mutations.txt")
  writeLines(log_no_mimp_res, con = "samples_without_kinase_rewiring_results.txt")
  
  return(full_results)
}  

run_mimp = function(fasta_path, phospho_path, search_engine, protein_id_col,
                    mut_data_path, ids_path, mutation_type_col, transcript_id_col, 
                    phosphosite_col, mutation_AA_change_colname){

  dir.create("mimp_results_dir")
  setwd("mimp_results_dir")
  dir.create("results_by_sample")

  # Loading NP data, different for different datatypes
  load(ids_path)
  
  # prepare fasta file for mimp input
  seqdata = format_fasta_file(fasta_path)

  # read phospho gct file, format rdesc and mat
  phospho_gct = parse.gctx(phospho_path)
  phos_cid = phospho_gct@cid %>% sub('^X', '', .)
  phos_rdesc = format_phospho_rdesc(phospho_gct, search_engine, protein_id_col, protein_id_type)
  phos_mat = format_phospho_mat(phospho_gct, phos_rdesc)
  
  # prepare mutation maf file before for loop sample-wise processing
  mut_data = format_mutation_full(mut_data_path, ids, mutation_type_col, transcript_id_col, seqdata)

  # prepare sample-wise mutation and phospho dfs, and then run mimp on each sample
  heatmap_df = run_mimp_samplewise(phos_cid, mut_data, seqdata, phos_rdesc, phos_mat,
                                   phosphosite_col, mutation_AA_change_colname)

  # create heatmap of predicted kinases results
  if(nrow(heatmap_df)>0){
    generate_mimp_heatmap(heatmap_df, groups_file_path, groups_file_SampleID_column)
  }
  
  return(heatmap_df)
}

mimp_results = run_mimp(fasta_path = fasta_path, 
                        phospho_path = phospho_path,
                        search_engine = search_engine, 
                        protein_id_col = protein_id_col,
                        mut_data_path = mut_data_path, 
                        ids_path = ids_path,
                        mutation_type_col = mutation_type_col, 
                        transcript_id_col = transcript_id_col,
                        phosphosite_col = phosphosite_col, 
                        mutation_AA_change_colname = mutation_AA_change_colname)