### Create Rmarkdown report for BlackSheep module ###

# tar_file  - URL of tar file created by task panoply_blacksheep
# yaml_file - URL of integrated yaml file produced by parameter_manager.r

args = commandArgs()

tar_file = args[1]
yaml_file = args[2]

library(yaml)
library(rmarkdown)
library(stringr)
library(dplyr)

rmd_blacksheep = function(tar_file, yaml_file){
  
  # extract values from final yaml file
  yaml_params = read_yaml(yaml_file)
  fdr_value = yaml_params$panoply_blacksheep$fdr_value
  SampleID_column = yaml_params$DEV_sample_annotation$sample_id_col_name
  groups_file = yaml_params$panoply_blacksheep$groups_file
  
  # extract files from tarball
  untar(tar_file)

  # begin writing rmd
  rmd = paste('---
title: "BlackSheep outlier analysis results"
output: 
  html_document:
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: true
      smooth_scroll: true
---\n
## Overview
This report summarizes the significant results (FDR <', fdr_value, ') of the BlackSheep module, which uses the blacksheepr package for differential extreme value analysis. Briefly, this module counts outliers in user-submitted data, tabulates outliers per group, and runs enrichment analysis (Fisher\'s exact test) to identify significant outliers in the groups of interest. More information can be found in Blumenberg et al. (2019, [preprint](https://www.biorxiv.org/content/10.1101/825067v2.full.pdf)).
              ')

  # if outlier analysis was performed (i.e. groups file provided), create outlier analysis rmd report
  if (length(grep("outlieranalysis", list.files("blacksheep", recursive = TRUE)))>0){
    for (i in c("positive", "negative")){
      rmd = paste0(rmd, '\n## ', str_to_title(i), ' outliers (FDR < ', fdr_value, ')')
      group_names = read.csv(groups_file) %>%
        colnames(.)
      group_names = group_names[!(group_names %in% SampleID_column)]
      for (group_name in group_names){
        rmd = paste(rmd, '\n###', group_name)
        outlier_filenames = grep(paste0("outlieranalysis_for_", group_name), list.files(file.path("blacksheep",i)), value=TRUE)
        for (file in outlier_filenames){
          outlier_analysis = read.csv(file.path("blacksheep", i, file))
          category = gsub(paste0(".*", group_name, "_(.+?)_.*"), "\\1", file)
          
          fdrcols = grep("fdr", colnames(outlier_analysis), value = TRUE)
          fdr_col = fdrcols[which(str_detect(fdrcols, "not", negate = TRUE))]
          
          outlier_analysis = outlier_analysis %>%
            select(gene, all_of(fdr_col)) %>%
            arrange_at(fdr_col) %>%
            rename(fdr = all_of(fdr_col)) %>%
            mutate(fdr = signif(fdr, 3))
          
          outlier_analysis = outlier_analysis[which(outlier_analysis[,"fdr"]< fdr_value),]
          
          save(outlier_analysis, file = paste(i, group_name, category, "outlier.RData", sep = "_"))
  
          rmd = paste0(rmd, '\n#### **', group_name, ': ', category, '**
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("', paste(i, group_name, category, "outlier.RData", sep = "_"), '")
```
Number of significant genes: `r dim(outlier_analysis)[1]`
\n
                      ')
          
          # if any significant genes were identified, show table + heatmap
          if (dim(outlier_analysis)[1] > 0){
            rmd = paste0(rmd, '\n
\n**Table**: Significant genes ranked by FDR value.
```{r echo=FALSE, warning=FALSE, message=FALSE}
library(DT)
datatable(outlier_analysis, rownames = FALSE, width = "500px")
```

**Figure**: Heatmap depicting the fraction of outliers within each significant gene for each sample. Sample annotations are ordered by ', group_name, '. Note: row names may not be annotated if there are more than 100 significant genes.
![](', file.path("blacksheep", i, grep(paste0(group_name, "_", category, ".png"), list.files(file.path("blacksheep", i)), value = TRUE)), ')


                       ')
          }
        }
      }
    }
  } else {
    rmd = paste(rmd, '\n## No outlier analysis was performed\n
No groups file was provided so enrichment analysis of outliers was not performed, only outlier count tables were calculated. See output tar file: blacksheep_outlier_analysis.tar.gz')
  }
  
  # write .rmd file
  writeLines(rmd, con = "rmd_blacksheep.rmd")
  
  # render .rmd file
  rmarkdown::render("rmd_blacksheep.rmd")
}

# run rmd_blacksheep function to make rmd report
rmd_blacksheep(tar_file, yaml_file)
