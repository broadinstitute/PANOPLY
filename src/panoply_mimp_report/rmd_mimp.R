### Create Rmarkdown report for mimp module ###

# tar_file  - URL of tar file created by task panoply_mimp
# output_prefix - prefix to be used for naming html file

# args = commandArgs(TRUE)
# 
# tar_file = as.character(args[1])
# output_prefix = as.character(args[2])

library(yaml)
library(rmarkdown)
library(stringr)
library(dplyr)

rmd_blacksheep = function(tar_file, output_prefix, type){
  
  # # extract files from tarball
  # untar(tar_file)
  
  # # extract values from final yaml file
  # yaml_params = read_yaml("final_output_params.yaml")
  
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
This report summarizes the results of the MIMP module, which uses the [rmimp R package](https://github.com/omarwagih/rmimp) for predicting the impact of mutations on kinase-substrate phosphorylation. Briefly, this module... More information about the MIMP algorithm can be found in [Wagih et al. 2015](https://www.nature.com/articles/nmeth.3396), and detailed documentation of the panoply_mimp module can be found here (insert link when WIKI page ready).
              ')
  
  
  
}