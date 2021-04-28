### Create Rmarkdown report for mimp module ###

# tar_file  - URL of tar file created by task panoply_mimp
# output_prefix - prefix to be used for naming html file

args = commandArgs(TRUE)

tar_file = as.character(args[1])
output_prefix = as.character(args[2])

library(yaml)
library(rmarkdown)
library(stringr)
library(dplyr)

rmd_mimp = function(tar_file, 
                    output_prefix){
  
  # extract files from tarball
  untar(tar_file)
  
  # begin writing rmd
  rmd = paste0('---
title: "MIMP: prediction of mutation impact on kinase-substrate phosphorylation for ', output_prefix, '"
output: 
  html_document:
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: true
      smooth_scroll: true
---\n
## Overview
This report summarizes the results of the panoply_mimp module, 
which uses the [rmimp R package](https://github.com/omarwagih/rmimp) for predicting 
the impact of missense mutations on kinase-substrate phosphorylation. 
Briefly, this module takes in mutation, phosphoproteomic, and protein sequence data 
and uses kinase specificity models to first identify mutations in close proximity to 
phosphosites (called phosphorylation-related SNV, or "pSNV") and then to predict if 
each pSNV disrupts an existing phosphosite or creates a new phosphosite, thereby 
potentially altering known kinase networks (referred to as a *kinase rewiring event*). More information about the MIMP algorithm 
can be found in [Wagih et al. 2015](https://www.nature.com/articles/nmeth.3396), 
and detailed documentation of the panoply_mimp module can be found here (insert link when WIKI page ready).

## MIMP results: predicted kinase rewiring
              ')
  
  if (file.exists("results_dir/mimp_output_kinase_rewiring_events_all.csv")){
    kinase_list = read.csv("results_dir/mimp_output_kinase_rewiring_events_all.csv") %>%
      mutate(score_wt = signif(as.numeric(score_wt), 3)) %>%
      mutate(score_mt = signif(as.numeric(score_mt), 3)) %>%
      mutate(log_ratio = signif(as.numeric(log_ratio), 3)) %>%
      mutate(prob = signif(as.numeric(prob), 3)) %>%
      select(-nseqs) %>%
      rename(kinase = kinase_pwm, kinase_fam = kinase_pwm_fam, 
             wt_psite_flanks = wild_type_seq, mt_psite_flanks = mut_seq)
    
    save(kinase_list, file = "kinase_list_edited.Rdata")
    
    rmd = paste0(rmd, '\n
### Full MIMP results

**Table**: list of all mutations predicted by MIMP to cause kinase rewiring events in each patient. Note that some mutations may affect multiple kinases, and some kinase activity may be affected by multiple mutations.
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("kinase_list_edited.Rdata")
library(DT)
datatable(kinase_list, rownames = FALSE)
```
*Column name key*:

* wt_psite_flanks = wild-type 7 amino acid flanking sequence of the phosphosite located at phosphosite_position
* mt_psite_flanks = mutated 7 amino acid flanking sequence of the phosphosite located at phosphosite_position
* mut_dist = distance of mutation from phosphosite central STY residue
* score_wt = a measure of how well a phosphosite with wild-type flanking sequence is likely phosphorylated by the given kinase (0 = no phosphorylation, 1 = perfect phosphorylation)
* score_mt = a measure of how well a phosphosite with a neighboring mutation is likely phosphorylated by the given kinase (0 = no phosphorylation, 1 = perfect phosphorylation)
* log_ratio = log(score_mt/score_wt) for the given kinase; 
positive log ratio indicates gain-of-phosphorylation event due to mutation and 
negative log ratio indicates loss-of-phosphorylation event, with magnitude indicating the confidence.
* prob = probability/confidence of the predicted kinase rewiring event
* effect = direction of predicted kinase rewiring: either "gain" or "loss" of phosphorylation by the given kinase

For more detailed description of the MIMP output values, please see the [MIMP help page](http://mimp.baderlab.org/help).


**Figure**: Heatmap depicting all predicted kinase rewiring events due to missense mutations. 
Each row represents a specific point mutation in a given gene/protein id and the predicted kinase it affects.
Heatmap values correspond to the log ratio of mutation/wild-type scores for a given kinase+gene+mutation combination.
Note: row names may not be annotated if there are more than 100 altered kinases predicted. See complete values underlying this heatmap at results_dir/kinase_rewiring_events_matrix_kinase_gene_mut_level.csv.
                 
![](results_dir/mimp_results_predicted_kinase_rewiring_kinase_gene_mut_level_heatmap.png)


\n
### Overall altered kinase activity due to missense mutations

```{r echo=FALSE, warning=FALSE, message=FALSE}
kinases = read.csv("results_dir/kinase_rewiring_events_matrix_kinase_level.csv", row.names = 1)
```
**Summary**: `r dim(kinases)[1]` kinases with altered activity due to mutations across `r dim(kinases)[2]` samples.


**Figure**: Heatmap depicting overall altered kinase activity in samples with predicted kinase rewiring events due to missense mutations. 
Heatmap values correspond to the log ratio of mutation/wild-type scores for a given kinase.
If a given kinase\'s activity was predicted to be affected by several mutations, that kinase is represented by its maximum absolute value log ratio. 
Note: row names may not be annotated if there are more than 100 altered kinases predicted. See complete values underlying this heatmap at results_dir/kinase_rewiring_events_matrix_kinase_level.csv.
                 
![](results_dir/mimp_results_predicted_kinase_rewiring_kinase_level_heatmap.png)
                ')


  } else {
    rmd = paste0(rmd, '\n
No predicted kinase rewiring events in any samples.')
  }
  
  rmd = paste0(rmd, '\n
## Phosphorylation-related SNVs\n

A phosphorylation-related SNV (pSNV) is defined as an SNV that occurs within +/- 7 amino acids of a phosphorylation site (determined from the input phosphoproteomic data). \n

               ')
  
  pSNV_tab = data.frame(patient_id = as.character(),
                        num_mut = as.numeric(),
                        num_pSNV = as.numeric(),
                        num_pSNV_remove_psite = as.numeric(),
                        num_pSNV_kinase_rewiring = as.numeric())
  muts_df = data.frame(patient_id = as.character(),
                      num_muts = as.numeric())
  for (sample in list.dirs("results_dir/results_by_sample/", recursive = FALSE, full.names = FALSE)){
    if (file.exists(file.path("results_dir", "results_by_sample", sample, "mimp_input", paste0(sample, "_mutation_mimp_input.csv")))){
      mut = read.csv(file.path("results_dir", "results_by_sample", sample, "mimp_input", paste0(sample, "_mutation_mimp_input.csv")))
      mut_sample = data.frame(patient_id = sample,
                              num_muts = as.numeric(nrow(mut)))
      muts_df = rbind(muts_df, mut_sample)
    }

    if (length(which(grepl("kinase_rewiring", list.files(file.path("results_dir", "results_by_sample", sample, "mutation_info"))))) == 1){
      pSNV = read.csv(grep("kinase_rewiring", list.files(file.path("results_dir", "results_by_sample", sample, "mutation_info"), full.names = T), value = T)) 
      pSNV_sample = data.frame(patient_id = sample,
                               num_mut = as.numeric(nrow(mut)),
                               num_pSNV = as.numeric(nrow(pSNV)),
                               num_pSNV_remove_psite = as.numeric(length(which(pSNV$mut_distance_from_phosphosite == 0))),
                               num_pSNV_kinase_rewiring = as.numeric(length(which(!is.na(pSNV$kinase_rewiring_event)))))
      pSNV_tab = rbind(pSNV_tab, pSNV_sample)
      
    } else if (length(which(grepl("pSNV_mutations", list.files(file.path("results_dir", "results_by_sample", sample, "mutation_info"))))) == 1) {
      pSNV = read.csv(grep("pSNV", list.files(file.path("results_dir", "results_by_sample", sample, "mutation_info"), full.names = T), value = T))
      pSNV_sample = data.frame(patient_id = sample,
                               num_mut = as.numeric(nrow(mut)),
                               num_pSNV = as.numeric(nrow(pSNV)),
                               num_pSNV_remove_psite = as.numeric(length(which(pSNV$mut_distance_from_phosphosite == 0))),
                               num_pSNV_kinase_rewiring = 0)
      pSNV_tab = rbind(pSNV_tab, pSNV_sample)

    }
    
  }
   
  save(muts_df, file = "all_mutations.Rdata")
  if (dim(pSNV_tab)[1]>0){
    save(pSNV_tab, file = "pSNV_table.Rdata")
    
    rmd = paste0(rmd, '\n
\n

```{r echo=FALSE, warning=FALSE, message=FALSE}
load("pSNV_table.Rdata")
load("all_mutations.Rdata")
```
**Summary**: \n

* Number of samples with at least one pSNV detected: `r dim(pSNV_tab)[1]` of `r dim(muts_df)[1]` total samples with mutations
* Number of samples with at least one pSNV that removes a central phosphosite STY residue: `r length(which(pSNV_tab$num_pSNV_remove_psite > 0))` of `r dim(muts_df)[1]` total samples with mutations 
* Number of samples with at least one pSNV predicted to cause a kinase rewiring event: `r length(which(pSNV_tab$num_pSNV_kinase_rewiring > 0))` of `r dim(muts_df)[1]` total samples with mutations \n 


**Table**: for each sample with at least one pSNV, breakdown of the total number of mutations in the sample, 
number of pSNVs, number of pSNVs affecting a central phosphosite residue, and 
number of pSNVs predicted by MIMP to rewire kinase activity.

```{r echo=FALSE, warning=FALSE, message=FALSE}
library(DT)
datatable(pSNV_tab, rownames = FALSE)
```

*Column name key*:

* patient_id = sample identifier 
* num_mut = total number of mutations
* num_pSNV = number of pSNVs
* num_pSNV_remove_psite = number of pSNVs that remove a central phosphosite STY residue, causing loss of a phosphosite
* num_pSNV_kinase_rewiring = number of pSNVs predicted by MIMP to cause kinase rewiring event(s)

For a complete list of pSNVs within each sample, please see the
pSNV_mutations_within_7AA_of_phosphosite.csv files located at 
results_dir/results_by_sample/[*sample name*]/mutation_info in the MIMP output tarball.
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
      left_join(muts_df) %>%
      mutate(percent_mismatch = signif(num_mismatch/num_muts,3)*100)
    
    save(mismatches, file = "mismatches_table.Rdata")
    
    rmd = paste0(rmd, '\n
##  Amino acid mismatches between mutation reference and fasta sequence input


*Please note that mutations mismatched with the input fasta sequence are removed from 
MIMP analysis. For a full table showing mismatches in each sample, 
please see results_dir/all_mismatchedAA_mutation_vs_fasta.csv.*

```{r echo=FALSE, warning=FALSE, message=FALSE}
load("mismatches_table.Rdata")
```

**Summary**:\n
* Number of samples with amino acid mismatches: `r dim(mismatches)[1]` of `r dim(muts_df)[1]` total samples with mutations
* Average percent of mismatched mutations per sample: `r signif(mean(mismatches$percent_mismatch, na.rm = TRUE),3)` %

**Figure**: for samples with amino acid mismatches between mutation reference and fasta sequence input, 
boxplot showing distribution of percent mismatched mutations.
```{r echo=FALSE, warning=FALSE, message=FALSE}
boxplot(mismatches$percent_mismatch, main = "Percent of mutations that are mismatched", ylab = "Percent")
```

**Table**: for samples with amino acid mismatches between mutation reference and fasta sequence input, 
breakdown of number of mismatches, number of total mutations, and percent of mutations that are mismatched.
```{r echo=FALSE, warning=FALSE, message=FALSE}
library(DT)
datatable(mismatches, rownames = FALSE)
```

*Column name key*:

* patient_id = sample identifier
* num_mismatch = number of mutations with mismatches between mutation reference and fasta sequence
* num_muts = total number of mutations
* percent_mismatch = percent of total mutations that are mismatched (unit = %)

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

rmd_mimp(tar_file,
  output_prefix)
