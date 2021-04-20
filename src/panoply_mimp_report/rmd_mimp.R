### Create Rmarkdown report for mimp module ###

# tar_file  - URL of tar file created by task panoply_mimp
# output_prefix - prefix to be used for naming html file

# args = commandArgs(TRUE)
# 
# tar_file = as.character(args[1])
# output_prefix = as.character(args[2])

setwd("C:/Users/karen/mimp/results_dir_brca/")
output_prefix = "brca"

library(yaml)
library(rmarkdown)
library(stringr)
library(dplyr)

rmd_mimp = function(tar_file, output_prefix){
  
  # # extract files from tarball
  # untar(tar_file)
  
  # begin writing rmd
  rmd = paste0('---
title: "Prediction of mutation impact on kinase-substrate phosphorylation using MIMP"
output: 
  html_document:
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: true
      smooth_scroll: true
---\n
## Overview
This report summarizes the results of the MIMP module, which uses the [rmimp R package](https://github.com/omarwagih/rmimp) for predicting the impact of missense mutations on kinase-substrate phosphorylation. Briefly, this module takes in mutation, phosphoproteomic, and protein sequence data and uses kinase specificity models to evaluate whether a mutation causes a phosphorylation-related SNV (pSNV) and to predict if the resulting pSNV disrupts an existing phosphosite or creates a new phosphosite, thereby potentially altering kinase networks. More information about the MIMP algorithm can be found in [Wagih et al. 2015](https://www.nature.com/articles/nmeth.3396), and detailed documentation of the panoply_mimp module can be found here (insert link when WIKI page ready).

## Overall altered kinase activity due to missense mutations
              ')
  
  if (file.exists("results_dir/mimp_output_kinase_rewiring_events_all.csv")){
    rmd = paste0(rmd, '\n
```{r echo=FALSE, warning=FALSE, message=FALSE}
kinases = read.csv("results_dir/kinase_rewiring_events_matrix_kinase_level.csv", row.names = 1)
```
**Summary**: `r dim(kinases)[1]` kinases with altered activity due to mutations across `r dim(kinases)[2]` samples.


**Figure**: Heatmap depicting overall altered kinase activity in samples with predicted kinase rewiring due to missense mutations close to a phosphosite. 
Heatmap values correspond to the log ratio of mutation vs wild type scores 
(i.e. how well a phosphosite is likely to be phosphorylated by the given kinase); 
positive log ratio indicates gain-of-phosphorylation event due to mutation and 
negative log ratio indicates loss-of-phosphorylation event, 
with magnitude indicating the confidence. 
If several mutations gave rise to altered activity in a given kinase, 
that kinase is represented by the maximum absolute value log ratio. 
Note: row names may not be annotated if there are more than 100 altered kinases predicted.
![](results_dir/mimp_results_predicted_kinase_rewiring_kinase_level_heatmap.png)

**Figure**: Heatmap depicting all predicted kinase rewiring events due to missense mutations close to a phosphosite. 
Each row represents a specific point mutation in a given gene/protein id and the predicted kinase acted on.
Heatmap values correspond to the log ratio of mutation vs wild type scores 
(i.e. how well a phosphosite is likely to be phosphorylated by the given kinase); 
positive log ratio indicates gain-of-phosphorylation event due to mutation and 
negative log ratio indicates loss-of-phosphorylation event, 
with magnitude indicating the confidence. 
Note: row names may not be annotated if there are more than 100 altered kinases predicted.
![](results_dir/mimp_results_predicted_kinase_rewiring_kinase_gene_mut_level_heatmap.png)
                                  ')
    
    kinase_list = read.csv("results_dir/mimp_output_kinase_rewiring_events_all.csv") %>%
      mutate(score_wt = signif(as.numeric(score_wt), 3)) %>%
      mutate(score_mt = signif(as.numeric(score_mt), 3)) %>%
      mutate(log_ratio = signif(as.numeric(log_ratio), 3)) %>%
      mutate(prob = signif(as.numeric(prob), 3)) %>%
      select(-c(mt, nseqs)) %>%
      rename(kinase = pwm, kinase_fam = pwm_fam, wt_psite_flanks = wt)
    
    save(kinase_list, file = "kinase_list_edited.Rdata")
    
    rmd = paste0(rmd, '\n
**Table**: full MIMP results showing all predicted kinase rewiring events per patient
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("kinase_list_edited.Rdata")
library(DT)
datatable(kinase_list, rownames = FALSE)
```
                 ')

  } else {
    rmd = paste0(rmd, 'No predicted kinase rewiring events in any samples.')
  }
  
  rmd = paste0(rmd, '\n
## Mutation information
               ')
  
  pSNV_tab = data.frame(patient_id = as.character(),
                        num_mut = as.numeric(),
                        num_pSNV = as.numeric(),
                        num_pSNV_remove_psite = as.numeric(),
                        num_pSNV_kinase_rewiring = as.numeric())
  muts_df = data.frame(patient_id = as.character(),
                      num_muts = as.numeric())
  for (sample in list.dirs("results_dir/results_by_sample/", recursive = FALSE)){
    mut = read.csv(file.path("results_dir", "results_by_sample", sample, "mimp_input", paste0(sample, "_mutation_mimp_input.csv")))
    mut_sample = data.frame(patient_id = sample,
                            num_muts = as.numeric(nrow(mut)))
    muts_df = rbind(muts_df, muts_sample)

    if (length(which(grepl("kinase_rewiring", list.files(file.path(sample, "mutation_info"))))) == 1){
      pSNV = read.csv(grep("kinase_rewiring", list.files(file.path(sample, "mutation_info"), full.names = T), value = T)) 
      pSNV_sample = data.frame(patient_id = sample,
                               num_mut = as.numeric(nrow(mut)),
                               num_pSNV = as.numeric(nrow(pSNV)),
                               num_pSNV_remove_psite = as.numeric(length(which(pSNV$mut_distance_from_phosphosite == 0))),
                               num_pSNV_kinase_rewiring = as.numeric(length(which(!is.na(pSNV$kinase_rewiring_event)))))
      
    } else if (length(which(grepl("pSNV_mutations", list.files(file.path(sample, "mutation_info"))))) == 1) {
      pSNV = read.csv(grep("pSNV", list.files(file.path(sample, "mutation_info"), full.names = T), value = T))
      pSNV_sample = data.frame(patient_id = sample,
                               num_mut = as.numeric(nrow(mut)),
                               num_pSNV = as.numeric(nrow(pSNV)),
                               num_pSNV_remove_psite = as.numeric(length(which(pSNV$mut_distance_from_phosphosite == 0))),
                               num_pSNV_kinase_rewiring = 0)

    }
    pSNV_tab = rbind(pSNV_tab, pSNV_sample)
  }
   
  save(muts_df, file = "all_mutations.Rdata")
  if (dim(pSNV_tab)[1]>0){
    save(pSNV_tab, file = "pSNV_table.Rdata")
    
    rmd = paste0(rmd, '\n
### Mutations in proximity to phosphosite and kinase rewiring

A phosphorylation-related SNV (pSNV) is defined as an SNV that occurs within +/- 7 amino acids of a phosphorylation site detected in the input data. \n

\n
**Summary**: \n

```{r echo=FALSE, warning=FALSE, message=FALSE}
load("pSNV_table.Rdata", "all_mutations.Rdata")
```

Number of samples with at least one pSNV detected: `r dim(pSNV_tab)[1]` of `r dim(muts_df)[1]` total samples \n 

Number of samples with at least one pSNV that removes a central phosphosite residue: `r length(which(pSNV_tab$num_pSNV_remove_psite > 0))` of `r dim(muts_df)[1]` total samples \n 

Number of samples with at least one pSNV predicted to cause kinase rewiring: `r length(which(pSNV_tab$num_pSNV_kinase_rewiring > 0))` of `r dim(muts_df)[1]` total samples \n 


**Table**: for each sample with at least one phosphorylation-related SNV (pSNV), breakdown of number of total mutations, number of pSNVs, and number of pSNVs predicted to alter kinase activity.
Column name key: patient_id = sample identifier;
num_mut = total number of mutations in the sample;
num_pSNV = number of pSNVs detected in the sample;
num_pSNV_remove_psite = number of pSNVs that cause loss of a phosphosite through change of central residue,
num_pSNV_kinase_rewiring = number of pSNVs predicted by MIMP to rewire kinase activity

```{r echo=FALSE, warning=FALSE, message=FALSE}
library(DT)
datatable(pSNV_tab, rownames = FALSE)
```

For a complete list of pSNVs within each sample and their effects, please see the pSNV_mutations_within_7AA_of_phosphosite.csv files within the mutation_info directories of each sample subdirectory.
                 ')
  } else {
    rmd = paste0(rmd, '\n
No phosphorylation-related SNVs detected in any samples.
                 ')
  }
  
  if (file.exists("results_dir/all_mismatchedAA_mutation_vs_fasta.csv")){
    
    mismatches = read.csv("results_dir/all_mismatchedAA_mutation_vs_fasta.csv") %>%
      group_by(patient_id) %>% 
      count(name = "num_mismatch") %>%
      mutate(num_mismatch = as.numeric(num_mismatch)) %>%
      left_join(muts) %>%
      mutate(percent_mismatch = num_mismatch/num_muts)
    
    save(mismatches, file = "mismatches_table.Rdata")
    
    rmd = paste0(rmd, '\n
###  Amino acid mismatches between mutation reference and fasta sequence input


*Please note that mutations mismatched with the input fasta sequence are removed from MIMP analysis. For a full table showing what the mismatches are, please see results_dir/all_mismatchedAA_mutation_vs_fasta.csv.*

```{r echo=FALSE, warning=FALSE, message=FALSE}
load("mismatches_table.Rdata")
```

Number of samples with amino acid mismatches: `r dim(mismatches)[1]` \n

Average percent of mismatched mutations per sample: `r mean(mismatches$percent_mismatch)` %
```{r echo=FALSE, warning=FALSE, message=FALSE}
boxplot(mismatches$percent_mismatch)
```

**Table**: for samples with amino acid mismatches between mutation reference and fasta sequence, breakdown of number of mismatches, number of total mutations, and percent of mutations that are mismatched
```{r echo=FALSE, warning=FALSE, message=FALSE}
library(DT)
datatable(mismatches, rownames = FALSE)
```
                 ')
  } else {
    rmd = paste0(rmd, '\n
No amino acid mismatches between mutation reference and fasta sequence input.
                 ')
  }

  
  rmd_name = paste(output_prefix, "mimp_rmd.rmd", sep="_")
  
  # write .rmd file
  writeLines(rmd, con = rmd_name)
  
  # render .rmd file
  rmarkdown::render(rmd_name)
  
}