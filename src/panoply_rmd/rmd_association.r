### Create Rmarkdown report for association module ###

# tar_file  - URL of tar file created by task panoply_association
# yaml_file - URL of integrated yaml file produced by parameter_manager.r
# label     - character, name of desired folder in tarball
# type      - character, data type

# args = commandArgs()
# 
# tar_file = args[1]
# yaml_file = args[2]
# label = args[3]
# type = args[4]

library(rmarkdown)
library(yaml)
library(cmapR)
library(tibble)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(plotly)
library(seriation)
library(glue)

tar_file = "C:/Users/karen/PhosphoDIA/Github/PANOPLY/src/panoply_rmd/acetylome-28a7ada4-4d04-46dd-8fe5-1c4e82449081-pgdac_main_full.tar"
yaml_file = "G:/Shared drives/Proteomics_LabMembers/LabMembers/Karen/blacksheep_scripts/test.yaml"
label = "panoply-lscc-3-2-acetylome-full-results"
type = "acetylome"

setwd("C:/Users/karen/PhosphoDIA/Github/PANOPLY/src/panoply_rmd")
source("heatmap_function.R")

hallmark_process_category <- c(
  HALLMARK_TNFA_SIGNALING_VIA_NFKB='signaling',
  HALLMARK_HYPOXIA='pathway',
  HALLMARK_CHOLESTEROL_HOMEOSTASIS='metabolic',
  HALLMARK_MITOTIC_SPINDLE='proliferation',
  HALLMARK_WNT_BETA_CATENIN_SIGNALING='signaling',
  HALLMARK_TGF_BETA_SIGNALING='signaling',
  HALLMARK_IL6_JAK_STAT3_SIGNALING='immune',
  HALLMARK_DNA_REPAIR='DNA damage',
  HALLMARK_G2M_CHECKPOINT='proliferation',
  HALLMARK_APOPTOSIS='pathway',
  HALLMARK_NOTCH_SIGNALING='signaling',
  HALLMARK_ADIPOGENESIS='development',
  HALLMARK_ESTROGEN_RESPONSE_EARLY='signaling',
  HALLMARK_ESTROGEN_RESPONSE_LATE='signaling',
  HALLMARK_ANDROGEN_RESPONSE='signaling',
  HALLMARK_MYOGENESIS='development',
  HALLMARK_PROTEIN_SECRETION='pathway',
  HALLMARK_INTERFERON_ALPHA_RESPONSE='immune',
  HALLMARK_INTERFERON_GAMMA_RESPONSE='immune',
  HALLMARK_APICAL_JUNCTION='cellular component',
  HALLMARK_APICAL_SURFACE='cellular component',
  HALLMARK_HEDGEHOG_SIGNALING='signaling',
  HALLMARK_COMPLEMENT='immune',
  HALLMARK_UNFOLDED_PROTEIN_RESPONSE='pathway',
  HALLMARK_PI3K_AKT_MTOR_SIGNALING='signaling',
  HALLMARK_MTORC1_SIGNALING='signaling',
  HALLMARK_E2F_TARGETS='proliferation',
  HALLMARK_MYC_TARGETS_V1='proliferation',
  HALLMARK_MYC_TARGETS_V2='proliferation',
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION='development',
  HALLMARK_INFLAMMATORY_RESPONSE='immune',
  HALLMARK_XENOBIOTIC_METABOLISM='metabolic',
  HALLMARK_FATTY_ACID_METABOLISM='metabolic',
  HALLMARK_OXIDATIVE_PHOSPHORYLATION='metabolic',
  HALLMARK_GLYCOLYSIS='metabolic',
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY='pathway',
  HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY='pathway',
  HALLMARK_P53_PATHWAY='proliferation',
  HALLMARK_UV_RESPONSE_UP='DNA damage',
  HALLMARK_UV_RESPONSE_DN='DNA damage',
  HALLMARK_ANGIOGENESIS='development',
  HALLMARK_HEME_METABOLISM='metabolic',
  HALLMARK_COAGULATION='immune',
  HALLMARK_IL2_STAT5_SIGNALING='signaling',
  HALLMARK_BILE_ACID_METABOLISM='metabolic',
  HALLMARK_PEROXISOME='cellular component',
  HALLMARK_ALLOGRAFT_REJECTION='immune',
  HALLMARK_SPERMATOGENESIS='development',
  HALLMARK_KRAS_SIGNALING_UP='signaling',
  HALLMARK_KRAS_SIGNALING_DN='signaling',
  HALLMARK_PANCREAS_BETA_CELLS='development'
)
hallmark_category = as.data.frame(hallmark_process_category)
names(hallmark_category) = "category"
hallmark_category = hallmark_category %>%
  rownames_to_column("pathway")

rmd_association = function(#tar_file, 
                           yaml_file, label, type){
  
  # extract values from final yaml file
  yaml_params = read_yaml(yaml_file)
  fdr_value = 0.01 # need to add from updated yaml file
  
  # extract files from tarball
  # untar(tar_file)
  
  assoc_dir = "association"
  gsea_dir = "ssgsea_assoc"
  
  # begin writing rmd
  rmd = paste0('---
title: "Association analysis results"
output: 
  html_document:
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: true
      smooth_scroll: true
---\n
## Overview
This report summarizes the results of the Association module at FDR < ', fdr_value, '.
              ')
  
  if (dir.exists(file.path(label, gsea_dir))){
    rmd = paste(rmd, '\n## ssGSEA Association Results
Please note that all volcano plots are interactive; hover mouse over a given point to see which pathway it corresponds to.
                ')
    
    for (in_dir in list.dirs(file.path(label, gsea_dir), full.names = FALSE, recursive = FALSE)){
      category = gsub(paste0("^ssgsea-", type, "-"), "", in_dir)
      category = gsub("-contrast_gc$", "", category)
      category_filename = gsub("\\.", "_", category)
      category = gsub("\\.", " ", category)
      rmd = paste(rmd, '\n###', category)
      
      combined_file = grep("combined", list.files(file.path(label, gsea_dir, in_dir), full.names = TRUE), value = TRUE)
      file = parse.gctx(combined_file)
      if (length(file@cid) == 1){
        mat = data.frame(file@mat) 
        comparison = gsub("\\.", " ", colnames(mat))
        colnames(mat) = paste("NES", colnames(mat), sep = "_")
        mat = mat %>%
          rownames_to_column("pathway")
        rdesc = data.frame(file@rdesc) %>%
          rownames_to_column("pathway") %>%
          rename(fdr = all_of(contains("fdr"))) %>%
          select(pathway, fdr) %>%
          left_join(mat) %>%
          arrange(fdr)
        
        # add hallmark categories if hallmark database used
        if(sum(file@rid %in% hallmark_category$pathway) > 0){
          rdesc = left_join(rdesc, hallmark_category)
        }
        
        rdesc_volc = rdesc %>%
          mutate(group = ifelse(fdr<fdr_value, "Significant", "Not significant")) %>%
          mutate(pathway = gsub("^HALLMARK_", "", pathway))
        
        save(rdesc_volc, file = paste0(category_filename, ".RData"))
        
        nes_name = grep("NES", colnames(rdesc), value = TRUE)
        
        rmd = paste0(rmd, '\n
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("', category_filename, '.RData")
if ("category" %in% names(rdesc_volc)){
  p = ggplot(rdesc_volc, aes(x = ', nes_name, ', y = -log10(fdr), colour = group, text = pathway, group = category)) +
    geom_point() +
    geom_hline(yintercept = -log10(', fdr_value, '), linetype = "dashed")
  ggplotly(p, tooltip = c("text", "category"))
} else {
  p = ggplot(rdesc_volc, aes(x = ', nes_name, ', y = -log10(fdr), colour = group, text = pathway)) +
    geom_point() +
    geom_hline(yintercept = -log10(', fdr_value, '), linetype = "dashed")
  ggplotly(p, tooltip = "text")
}

```
**Figure**: Volcano plot summarizing ssGSEA pathway results for ', category, ', comparison ', comparison, '. X axis represents the Normalized Enrichment Score (NES) for ', comparison, '; positive NES values indicate enrichment in ', gsub(" over.*", "", comparison), ' and negative NES values indicate enrichment in ', gsub(".*over ", "", comparison), '. Y axis represents the -log10 of the FDR value. Results above the dashed line are significant at FDR cutoff = ', fdr_value, '.
                    ')
        
        rdesc_tab = rdesc %>% 
          filter(fdr < fdr_value) %>%
          mutate(fdr = signif(fdr, 3))
        
        if (dim(rdesc_tab)[1]>= 1){
          save(rdesc_tab, file = paste0(category_filename, "_filtered.RData"))
          rmd = paste0(rmd, '\n
**Table**: ', dim(rdesc_tab)[1], ' significantly enriched pathways with FDR < ', fdr_value, '.
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("', category_filename, '_filtered.RData")
library(DT)
datatable(rdesc_tab, rownames = FALSE, width = "500px")
```
                       ')
        } else if (dim(rdesc_tab)[1] == 0){
          rmd = paste0(rmd, '\n
No significantly enriched pathways with FDR < ', fdr_value, '.                       
                       ')
        }
      } else {
        mat = data.frame(file@mat)
        rdesc = file@rdesc
        
        rdesc2 = rdesc %>%
          select(contains("fdr")) %>%
          rownames_to_column("pathway") %>%
          gather(key = "comparison", value = "fdr", -pathway) %>%
          mutate(comparison = gsub("fdr\\.pvalue\\.", "", comparison)) %>%
          mutate(fdr = signif(fdr, 3))
        
        mat2 = mat %>%
          rownames_to_column("pathway") %>%
          gather(key = "comparison", value = "NES", -pathway, na.rm = TRUE) %>%
          left_join(rdesc2) %>%
          filter(fdr < fdr_value) %>%
          arrange(fdr)
        
        if(sum(rid %in% hallmark_category$pathway) > 0){
          mat2 = mat2 %>%
            left_join(hallmark_category) %>%
            arrange(category)
        }
   
        if (dim(mat2)[1]>=1){
          save(mat2, file = paste0(category_filename, "_filtered.RData"))
          
          pw_hm(output.prefix = file, fdr.max = fdr_value, n.max = 50, ptmsigdb=F)
          file.rename(paste0("heatmap_max.fdr_", fdr_value, "_n.max_50.png"), paste0(category_filename, "_heatmap_max.fdr_", fdr_value, "_n.max_50.png"))
          
          rmd = paste0(rmd, '\n
![Figure: heatmap](', category_filename, '_heatmap_max.fdr_', fdr_value, '_n.max_50.png)
\n
**Table**: ', dim(mat2)[1], ' significantly enriched pathways across ', length(file@cid), ' comparisons with FDR < ', fdr_value, '.
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("', category_filename, '_filtered.RData")
library(DT)
datatable(mat2, rownames = FALSE)
```
                       ')
        } else if (dim(mat2)[1] == 0) {
          rmd = paste0(rmd, '\n
No significantly enriched pathways across ', length(file@cid), ' comparisons with FDR < ', fdr_value, '.                       
                       ')
        }

      }
    }
  }
  
  # write .rmd file
  writeLines(rmd, con = "rmd_association.rmd")
  
  # render .rmd file
  rmarkdown::render("rmd_association.rmd")
  
}

# run rmd_association function to make rmd report
rmd_association(#tar_file, 
                yaml_file = yaml_file, label = label, type = type)
