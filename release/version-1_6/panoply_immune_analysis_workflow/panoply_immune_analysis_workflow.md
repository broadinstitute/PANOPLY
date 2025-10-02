# ```panoply_immune_analysis_workflow```

## Description

This workflow performs [immune analysis](./Data-Analysis-Modules%3A-panoply_immune_analysis) on transcriptomics data, assigning immune scores for understanding the tumor microenvironment; it also generates an [interactive report](./Report-Modules%3A-panoply_immune_analysis_report.md). It executes the following modules:

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_immune_analysis</code>](./Data-Analysis-Modules%3A-panoply_immune_analysis)         |  performs immune analysis  |
| [<code>panoply_immune_analysis_report</code>](./Report-Modules%3A-panoply_immune_analysis_report.md)   |  creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the immune analysis results |


## Input

### Required inputs:

* ```inputData```: (File) tarball output (`.tar`) for one of the other modules, or RNAseq data in `gct` format (when `standalone` is `TRUE`)
* ```type```: (String) proteome data type
* ```standalone```: (String) set to ```TRUE``` to run as a self-contained module; if ```TRUE``` the ```analysisDir``` input is required
* ```yaml```: (`.yaml` file) parameters in `yaml` format
* ```label```: (String) directory name within the tar file

### Optional inputs:

* ```analysisDir```: (String) name of analysis directory
* ```groupsFile```: (`.csv` file) subset of sample annotations, providing classes for enrichment analysis
* ```fdr```: (Float, default = 0.05) FDR cutoff value for significance in enrichment analysis
* ```heatmapWidth```: (Int, default = 10) width of heatmap generated 
* ```heatmapHeight```: (Int, default = 15) height of heatmap generated


## Output

`panoply_immune_analysis_workflow` produces the follow outputs:

* ```outputs```: (`.tar`) An output tar file that contains the immune analysis results 
* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
