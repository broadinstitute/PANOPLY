# ```panoply_mimp_report```

## Description

This module creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the results of the [panoply_mimp](./Data-Analysis-Modules%3A-panoply_mimp) module. An interactive table displays the full list of kinase rewiring events identified by MIMP in all samples, followed by two heatmaps: one depicting all kinase rewiring events predicted by MIMP due to missense mutations, and the other summarizing overall altered kinase activity in each sample as predicted by MIMP. Information about phosphorylation-related SNVs (missense mutations in close proximity to phosphosites) detected in the dataset is then summarized, followed by a report of any samples that have amino acid mismatches between the reference amino acid in the mutation data and the amino acid at this position in the fasta sequence input.

## Input

### Required inputs:

* ```tar_file```: (`.tar` file) output tar file from ```panoply_mimp```
* ```output_prefix```: (String) prefix used to name the output tar file in ```panoply_mimp```

## Output

* ```report_out```: (`.html` file) Interactive [R Markdown](https://rmarkdown.rstudio.com/) report.