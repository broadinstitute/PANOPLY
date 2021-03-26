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
# result_dir = as.character(args[2])
# yaml_file = as.character(args[3])

setwd("C:/Users/karen/mimp")
mut_data_path = "S:/CPTAC3/PGDAC/brca/prospective/v5.4-public/data-freeze/prosp-brca-v5.4-public-BRCA-freeze-v5.final_analysis_set.maf.txt"
phospho_path = "S:/CPTAC3/PGDAC/brca/prospective/v5.4-public/data-freeze/prosp-brca-v5.4-public-phosphoproteome-ratio-norm-NArm.gct"
ids_path = "P:/LabMembers/Abhijeet Mavi/MIMP_final/ids.RData"
fasta_path = "S:/CPTAC2/RefSeq_20160914/RefSeq.20160914_Human_ucsc_hg19_customProDBnr_mito_150contams.fasta"
report_type = "SpectrumMill"

## mutation file column names
protein_change_colname = "Protein_Change"
#protein_change_colname = "HGVSp_Short"
mutation_type_col = "Variant_Classification"
patient_id_col = "Tumor_Sample_Barcode"
RefSeq_col = "refseq_mrna_id"

## phospho file column names
report_type = "SpectrumMill"
if (report_type == "SpectrumMill"){
  phosphosite_col = "variableSites"
  accession_number_col = "accession_number"
} 
  

dir.create("results_dir")
setwd("results_dir")
dir.create("results_by_sample")


# Loading NP data, different for different datatypes
load(ids_path)


# Preprocessing the fasta file in RefSeq format
fasta = Biostrings::readAAStringSet(fasta_path)
names(fasta) <- names(fasta) %>% sub('\\..*', '', .) %>% sub('^>','', .)
keep = grepl("^NP", names(fasta))
seqdata = as.list( as.character( fasta[keep] ) )


# read phospho gct file and remove Xs in front of patient IDs
phospho_gct = parse.gctx(phospho_path)
phos_cdesc <- phospho_gct@cdesc
rownames(phos_cdesc) = sub('^X','', rownames(phos_cdesc))
phos_rdesc <- phospho_gct@rdesc
phos_mat <- phospho_gct@mat
colnames(phos_mat) <- sub('^X','', colnames(phos_mat))
phos_cid <- phospho_gct@cid %>% sub('^X', '', .)


# filter for fully localized sites (# localized = # actual sites)
if (report_type == "SpectrumMill"){
  phos_rdesc2 = phos_rdesc %>%
    filter(Best_numActualVMSites_sty == Best_numLocalizedVMsites_sty)
  rownames(phos_rdesc2) = phos_rdesc2$id
}

phos_mat2 = phos_mat[rownames(phos_rdesc2),]

phos_rdesc2$NP_revised = sub('\\..*','',phos_rdesc2[, accession_number_col])

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

log_no_mut = "Samples with no mutations (mimp not run): \n"
log_no_psnv = "Samples with no pSNV's: \n"
log_no_rewiring = "Samples with no kinase rewiring events: \n"
log_all_filtered = "Samples with no phosphorylation or sequence data remaining after filtering: \n"

full_results = data.frame(pwm = as.character(),
                          log_ratio = as.numeric(),
                          patient_id = as.character())

gain_site_all = data.frame(NP_id = as.character(),
                           mutation = as.character(),
                           patient_id = as.character())
loss_site_all = data.frame(NP_id = as.character(),
                           mutation = as.character(),
                           patient_id = as.character())

# Algorithm to run MIMP on whole dataset, patient-wise
for (i in phos_cid){
  
  dir.create(i)
  dir.create("mimp_results_per_sample")
  dir.create("mimp_input")
  dir.create("mimp_results")
  
  # select patient-specific mutation entries
  mut_i = mut_data3[which(grepl(i,mut_data3[, patient_id_col])),]
  write.csv(mut_i, file.path(i, "mimp_input", paste0(i, "_mutation_data_all.csv")), row.names = FALSE)
  
  if(dim(mut_i)[1] > 0){
    
    # Creating the dataframe for mutation file in MIMP format from respective columns of mutation file
    df_mut = mut_i %>%
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
    write.csv(df_mut, file.path(i, "mimp_input", paste0(i, "_mutation_mimp_input.csv")), row.names = FALSE)
    
    gain_site = df_mut2 %>%
      filter(!(str_sub(mutation, 1, 1) %in% c("S", "T", "Y"))) %>%
      filter(str_sub(mutation, -1) %in% c("S", "T", "Y"))
    if (dim(gain_site)[1] > 0){
      write.csv(gain_site, file.path(i, "mimp_results", paste0(i, "_central_residue_STY_gain.csv")), row.names = FALSE)
    }
    gain_site = gain_site %>% mutate(patient_id = i)
    gain_site_all = rbind(gain_site_all, gain_site)
    
    lose_site = df_mut2 %>%
      filter(!(str_sub(mutation, -1) %in% c("S", "T", "Y"))) %>%
      filter(str_sub(mutation, 1, 1) %in% c("S", "T", "Y"))
    if (dim(lose_site)[1] > 0){
      write.csv(lose_site, file.path(i, "mimp_results", paste0(i, "_central_residue_STY_loss.csv")), row.names = FALSE)
    }
    lose_site = lose_site %>% mutate(patient_id = i)
    loss_site_all = rbind(loss_site_all, lose_site)
    
    # save list of mutation reference AAs that don't match with fasta sequence
    # save list of mutations that fall in +/- 7 AA window of STY
    AA_mismatch = data.frame(protein = as.character(),
                             mutation = as.character(),
                             mutation_position = as.numeric(),
                             mutation_reference_AA = as.character(),
                             fasta_AA = as.character())
    AA_windows = data.frame(protein = as.character(),
                            mutation = as.character(),
                            window_containing_STY = as.character())
    for (j in 1:length(df_mut2$NP_id)){
      AA_pos = as.numeric(gsub("[A-Z]", "", df_mut2$mutation[j]))
      if (str_sub(seqdata[df_mut2$NP_id[j]], AA_pos, AA_pos) != str_sub(df_mut2$mutation[j], 1, 1)){
        j_df = data.frame(protein = df_mut2$NP_id[j],
                          mutation = df_mut2$mutation[j],
                          AA_position = AA_pos,
                          AA_mutation_reference = str_sub(df_mut2$mutation[j], 1, 1),
                          AA_fasta = str_sub(seqdata[df_mut2$NP_id[j]], AA_pos, AA_pos))
        AA_mismatch = rbind(AA_mismatch, j_df)
      }
      AA_window = str_sub(seqdata[df_mut2$NP_id[j]], AA_pos - 7, AA_pos + 7)
      if (grepl("S|T|Y", AA_window)){
        new_df = data.frame(protein = df_mut2$NP_id[j],
                            mutation = df_mut2$mutation[j],
                            window_containing_STY = AA_window)
        AA_windows = rbind(AA_windows, new_df)
      }
    }
    
    if (dim(AA_mismatch)[1] > 0){
      write.csv(AA_mismatch, file.path(i, "mimp_input", paste0(i, "_mismatchedAA_mutation_vs_fasta.csv")), row.names = FALSE)
    }
    
    if (dim(AA_windows)[1] > 0){
      write.csv(AA_windows, file.path(i, "mimp_results", paste0(i, "_mutations_with_neighboring_STY.csv")), row.names = FALSE)
    }
  
    # Creating the phospho data in MIMP format
    phos_mat_id = !is.na(phos_mat2[,i])
    
    df_phospho = phos_rdesc2[phos_mat_id,] %>%
      select(NP_revised,
             all_of(phosphosite_col))
    names(df_phospho)[2] = "phosphosite"
    df_phospho = df_phospho %>%
      mutate(phosphosite = gsub(' $','', phosphosite)) %>%
      mutate(phosphosite = strsplit(phosphosite, " ")) %>%
      unnest(phosphosite) %>%
      filter(grepl("S|T|Y|s|t|y", phosphosite)) %>%
      mutate(phosphosite = gsub("S|T|Y|s|t|y", "", phosphosite))
    write.csv(df_phospho, file.path(i, "mimp_input", paste0(i, "_phospho_mimp_input.csv")), row.names = FALSE)
    
    # note any protein identifiers from phospho data that are not in the fasta file
    prot_not_in_fasta = list(df_phospho$NP_revised[which(!(df_phospho$NP_revised %in% names(fasta)))])
    if (length(prot_not_in_fasta[[1]]) > 0){
      names(prot_not_in_fasta) = "phospho identifiers without fasta match"
      write.csv(prot_not_in_fasta, file.path(i, "mimp_input", paste0(i, "_phospho_identifiers_not_in_fasta.csv")), row.names = FALSE)
    }
    
    # Run MIMP (patient-specific)
    warning(paste0("Running MIMP for ", i))
    results_i = mimp(muts = df_mut2, seqs = seqdata, psites = df_phospho)
    
    if (!is.null(results_i)){
      write.csv(results_i, file = file.path(i, "mimp_results", paste(i, "mimp_output_kinase_rewiring_events.csv", sep = "_")), row.names = FALSE)
      
      res_filt = results_i %>%
        distinct() %>%
        select(pwm, log_ratio) %>%
        mutate(patient_id = i)
      
      full_results = rbind(full_results, res_filt)
        
    } else {
      if (length(which(grepl("No pSNVs were found!", names(last.warning))))>0){
        log_no_psnv = paste(log_no_psnv, i, sep = "\n")
      } else if (length(which(grepl("No phosphorylation or sequence data remaining after filtering!", names(last.warning))))>0){
        log_all_filtered = paste(log_all_filtered, i, sep = "\n")
      } else if (length(which(grepl("No rewiring events found, returning NULL", names(last.warning))))>0){
        log_no_rewiring = paste(log_no_rewiring, i, sep = "\n")
      }
    }
    
  } else {
    log_no_mut = paste(log_no_mut, i, sep = "\n")
  }
  
}

write.csv(full_results, "mimp_output_kinase_rewiring_events_all.csv", row.names = FALSE)
write.csv(gain_site_all, "all_central_residue_STY_gain.csv", row.names = FALSE)
write.csv(loss_site_all, "all_central_residue_STY_loss.csv", row.names = FALSE)
writeLines(log_no_mut, con = "samples_without_mutations.txt")
writeLines(log_no_psnv, con = "samples_without_pSNVs.txt")
writeLines(log_no_rewiring, con = "samples_without_kinase_rewiring.txt")
writeLines(log_all_filtered, con = "samples_all_data_filtered.txt")

full_results_edit = full_results %>%
  filter(!is.na(log_ratio)) %>%
  group_by(pwm, patient_id) %>%
  slice(which.max(abs(log_ratio))) %>%
  spread(key = pwm, value = log_ratio) %>%
  column_to_rownames("patient_id") %>%
  t()
results.max = max(full_results$log_ratio, na.rm = TRUE)
results.min = min(full_results$log_ratio, na.rm = TRUE)

color = colorRamp2(c(results.min, 0, results.max), c("blue", "white", "red"))

heatmap_all = Heatmap(full_results_edit, cluster_rows = FALSE, cluster_columns = FALSE, 
                      col = color, column_title = "Sample", row_title = "Kinase",
                      column_title_side = "bottom", column_title_gp = gpar(fontsize = 10),
                      row_title_side = "left", row_title_gp = gpar(fontsize = 10),
                      row_names_side = "left", row_names_gp = gpar(fontsize=6),
                      column_names_gp = gpar(fontsize = 6), 
                      heatmap_legend_param = list(title = "log_ratio"),
                      width = ncol(full_results_edit)*unit(2, "mm"),
                      height = nrow(full_results_edit)*unit(2, "mm"))

# print heatmap
pdf('mimp_results_heatmap.pdf')

draw(heatmap_all)
# barplot(t(unmatched_NP),xlab = 'Patient_ID', ylab = 'Percentage of RefSeq IDs matched')

dev.off()
