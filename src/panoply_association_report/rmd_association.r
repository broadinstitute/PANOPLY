### Create Rmarkdown report for association module ###

# tar_file  - URL of tar file created by task panoply_association
# yaml_file - URL of integrated yaml file produced by parameter_manager.r
# label     - character, name of desired folder in tarball
# type      - character, data type

args = commandArgs(TRUE)

tar_file = args[1]
yaml_file = args[2]
label = args[3]
type = args[4]

library(rmarkdown)
library(yaml)
library(cmapR)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(glue)

source("/prot/proteomics/Projects/PGDAC/src/rmd-ssgsea-functions.R")

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

volcano_plot = function(mat, rdesc, rid, split_word, fdr_value, category, category_filename, rmd){
  comparison = gsub("\\.", " ", colnames(mat))
  comparison_filename = gsub("\\.", "_", colnames(mat))
  colnames(mat) = paste("NES", colnames(mat), sep = "_")
  
  mat2 = mat %>%
    rownames_to_column("pathway")
  
  rdesc2 = rdesc %>%
    rownames_to_column("pathway") %>%
    rename(fdr = all_of(contains("fdr"))) %>%
    select(pathway, fdr) %>%
    left_join(mat2) %>%
    arrange(fdr)
  
  # add hallmark categories if hallmark database used
  if(sum(rid %in% hallmark_category$pathway) > 0){
    rdesc2 = left_join(rdesc2, hallmark_category)
  }
  
  rdesc_volc = rdesc2 %>%
    mutate(group = ifelse(fdr<fdr_value, "Significant", "Not significant")) %>%
    mutate(pathway = gsub("^HALLMARK_", "", pathway))
  
  save(rdesc_volc, file = paste0(category_filename, comparison_filename, ".RData"))
  
  nes_name = grep("NES", colnames(rdesc2), value = TRUE)
  
  if (split_word == "vs"){
    rmd = paste(rmd, '\n#### Contrast:', comparison)
  }
  
  rmd = paste0(rmd, '\n
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("', category_filename, comparison_filename, '.RData")
if ("category" %in% names(rdesc_volc)){
  p = ggplot(rdesc_volc, aes(x = ', nes_name, ', y = -log10(fdr), colour = group, text = pathway, group = category)) +
    geom_point() +
    geom_hline(yintercept = -log10(', fdr_value, '), linetype = "dashed") +
    scale_color_manual(values=c("black", "red"))
  ggplotly(p, tooltip = c("text", "category"))
} else {
  p = ggplot(rdesc_volc, aes(x = ', nes_name, ', y = -log10(fdr), colour = group, text = pathway)) +
    geom_point() +
    geom_hline(yintercept = -log10(', fdr_value, '), linetype = "dashed") +
    scale_color_manual(values=c("black", "red"))
  ggplotly(p, tooltip = "text")
}
```
**Figure**: Volcano plot summarizing ssGSEA pathway results for ', category, ', contrast ', comparison, '. X axis represents the Normalized Enrichment Score (NES) for ', comparison, '; positive NES values indicate enrichment in ', gsub(paste0(" ", split_word, ".*"), "", comparison), ' and negative NES values indicate enrichment in ', gsub(paste0(".*", split_word, " "), "", comparison), '. Y axis represents the -log10 of the FDR value. Results above the dashed line are significant at FDR cutoff = ', fdr_value, '.
               ')
        
  rdesc_tab = rdesc2 %>% 
    filter(fdr < fdr_value) %>%
    mutate(fdr = signif(fdr, 3))
  
  # if there are any significant results, make a table
  if (dim(rdesc_tab)[1]>= 1){
    save(rdesc_tab, file = paste0(category_filename, comparison_filename, "_filtered.RData"))
    rmd = paste0(rmd, '\n
**Table**: ', dim(rdesc_tab)[1], ' significantly enriched pathways for ', comparison, ' with FDR < ', fdr_value, '.
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("', category_filename, comparison_filename, '_filtered.RData")
library(DT)
datatable(rdesc_tab, rownames = FALSE, width = "500px")
```
\n
                 ')
  } else {
    rmd = paste0(rmd, '\n
No significantly enriched pathways for ', comparison, ' with FDR < ', fdr_value, '.  
\n
                 ')
  }
  
  return(rmd)
}

rmd_association = function(tar_file, yaml_file, label, type){
  
  # extract values from final yaml file
  yaml_params = read_yaml(yaml_file)
  fdr_value = yaml_params$panoply_association_report$fdr_value
  gsea_dir = "ssgsea_assoc"
  
  # extract files from tarball
  untar(tar_file)
  
  # begin writing rmd
  rmd = paste0('---
title: "Association analysis results"
output: 
  html_document:
    toc: true
    toc_depth: 4
    toc_float:
      collapsed: true
      smooth_scroll: true
---\n
## Overview
This report summarizes the results of the Association module at FDR < ', fdr_value, '. The association module runs several marker selection algorithms on each comparison pair to assign a combined score identifying the best markers. ssGSEA is then run on these marker results to identify significantly enriched pathways.
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
      
      # if there is only one comparison in GCT, make volcano plot
      if (length(file@cid) == 1){
        rmd = volcano_plot(mat = data.frame(file@mat), rdesc = data.frame(file@rdesc), rid = file@rid, split_word = "over", fdr_value = fdr_value, category = category, category_filename = category_filename, rmd = rmd)
        
      # if there are multiple comparisons in GCT, make volcano plots + heatmap  
      } else {
        mat = data.frame(file@mat)
        compare = colnames(mat)
        rdesc = data.frame(file@rdesc)
        
        for (name in colnames(mat)){
          mat_1 = mat %>%
            select(all_of(name))
          
          rdesc_1 = rdesc %>%
            select(contains(name))
          
          rmd = volcano_plot(mat = mat_1, rdesc = rdesc_1, rid = file@rid, split_word = "vs", fdr_value = fdr_value, category = category, category_filename = category_filename, rmd = rmd)
        }
      
        pw_hm(output.prefix = file, fdr.max = fdr_value, ptmsigdb=F)
        file.rename(paste0("heatmap_max.fdr_", fdr_value, "_n.max_10.png"), paste0(category_filename, "_heatmap_max.fdr_", fdr_value, "_n.max_10.png"))
        
        rmd = paste0(rmd, '\n#### Overview of all contrasts
![**Figure**: Heatmap summarizing significant ssGSEA pathway results for ', category, ', all contrasts, clustered by Hallmark process category. Asterisk denotes a significant result at FDR cutoff = ', fdr_value, '.](', category_filename, '_heatmap_max.fdr_', fdr_value, '_n.max_10.png)
\n                       
                     ')
      }
    }
  } else {
    rmd = paste(rmd, '\n## ssGSEA Association Results
No ssGSEA results were found.
                ')
  }

  rmd_name = paste(label, type, "association_rmd.rmd", sep = "_")
  
  # write .rmd file
  writeLines(rmd, con = rmd_name)
  
  # render .rmd file
  rmarkdown::render(rmd_name)
  
}

# run rmd_association function to make rmd report
rmd_association(tar_file, yaml_file = yaml_file, label = label, type = type)
