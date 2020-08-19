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

rmd_association = function(tar_file, 
                           yaml_file){
  
  # extract values from final yaml file
  yaml_params = read_yaml(yaml_file)
  fdr_value = 0.01 # need to add from updated yaml file
  
  # extract files from tarball
  untar(tar_file)
  
  assoc_dir = "association"
  gsea_dir = "ssgsea_assoc"
  
  # begin writing rmd
  rmd = paste('---
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
This report summarizes the significant results (FDR <', fdr_value, ') of the Association module....
              ')
  
  if (dir.exists(file.path(label, gsea_dir))){
    rmd = paste(rmd, '\n## ssGSEA Association Results')
    
    for (in_dir in list.dirs(file.path(label, gsea_dir), full.names = FALSE, recursive = FALSE)){
      category = gsub(paste0("^ssgsea-", type, "-"), "", in_dir)
      category = gsub("-contrast_gc$", "", category)
      category = gsub("\\.", " ", category)
      rmd = paste(rmd, '\n###', category)
      
      combined_file = grep("combined", list.files(file.path(label, gsea_dir, in_dir), full.names = TRUE), value = TRUE)
      file = parse.gctx(combined_file)
      if (length(file@cid) == 1){
        mat = data.frame(file@mat) 
        colnames(mat) = paste("NES", colnames(mat), sep = "_")
        mat = mat %>%
          rownames_to_column("rowname")
        rdesc = data.frame(file@rdesc) %>%
          rownames_to_column("rowname") %>%
          rename(fdr = all_of(contains("fdr"))) %>%
          select(rowname, fdr) %>%
          left_join(mat) %>%
          mutate(group = ifelse(fdr<fdr_value, "Significant", "Not significant"))
        
        save(rdesc, file = paste0(category, ".RData"))
        
        nes_name = grep("NES", colnames(rdesc), value = TRUE)
        
        rmd = paste0(rmd, '\n
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("', category, '.RData")
p = ggplot(rdesc, aes(x = ', nes_name, ', y = -log10(fdr), colour = group, text = rowname)) +
  geom_point() +
  geom_hline(yintercept = -log10(', fdr_value, '), linetype = "dashed")
ggplotly(p, tooltip = "text")
```                    
                    ')
      } else {
        mat = data.frame(file@mat)
        rdesc = data.frame(file@rdesc) %>%
          select(contains("fdr"))
        
        pw_hm(output.prefix = file, fdr.max = fdr_value, ptmsigdb=F)
        
        save(mat, rdesc, file = paste0(category, ".RData"))
        
        rmd = paste0(rmd, '\n
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("', category, '.RData")
morpheus(mat, rowAnnotations = rdesc, height = "1000px", na.rm = TRUE)
pheatmap(mat, height = "1000px")
``` 

\n\n\n
                     ')
        
      }
    }
  }
  
  # write .rmd file
  writeLines(rmd, con = "rmd_association.rmd")
  
  # render .rmd file
  rmarkdown::render("rmd_association.rmd")
  
}
