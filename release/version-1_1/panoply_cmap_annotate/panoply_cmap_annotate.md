# ```panoply_cmap_annotate```

## Description
Annotates CMAP output from `panoply_cmap_analysis` by performing the following analyses. This module is already integrated into the `panoply_cmap_analysis` module.

* To identify how many trans-correlated genes for all candidate regulatory genes can be directly explained by gene expression changes measured in the CMAP shRNA perturbation experiments, knockdown gene expression consensus signature z-scores (knockdown/control) are used to identify regulated genes, followed by counting the number of trans-genes in this list of regulated genes.

* To obtain biological insight into the list of candidate driver genes, enrichment analysis is performed on samples with extreme CNA values (amplification or deletion) to identify statistically enriched sample annotation subgroups from annotation provided in the `cmap_enrichment_groups` file.

* Prepares input data for subsequent GSEA analysis on cis/trans-correlation values to find enriched pathways.


## Input

Required inputs:

* ```tarball```: (`.tar` file) tarball from `panoply_cmap_connectivity` (internal task)
* ```cmap_data_file```: (String) CMAP level 5 gene knockdown data table in `gctx` format (gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/annotated_GSE92742_Broad_LINCS_Level5_COMPZ_geneKDsubset_n36720x12328.gctx)
* ```yaml```: (File) parameters file in `yaml` format


Optional inputs:

* ```cmap_grp```: (String, default = "all") filename prefix used for files created during CMAP analysis
* ```cmap_typ```: (String, default = "pome") omics data type; suppoted types are "pome" (all proteomics/PTM types) and "mrna" (for RNAseq data)
* ```cmap_enrichment_groups```: (`.csv` File) subset of sample annotations, providing classes for enrichment analysis of candidate genes
* ```cna_threshold```: (Float, default = 0.3) copy number up/down threshold; copy number is considered UP regulated if > `cna_threshold` and DOWN regulated if < `-cna_threshold` 
* ```log_transform```: (String, default = "FALSE") if TRUE, log2 transform input data
* ```alpha```: (String, default = 0.05) p-value threshold for CMAP profile z-scores (see Description) and enrichments
* ```outFile```: (String, default = "panoply_cmap-annotate-output.tar") output `.tar` file name

## Output
Tarball with annotated CMAP analysis results. See [Data Analysis Modules: panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis) for details.
