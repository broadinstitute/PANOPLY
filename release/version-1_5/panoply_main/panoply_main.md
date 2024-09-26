# ```panoply_main```

## Description
Performs proteogenomic analysis for a given ome (proteome, phosphoproteome, ubiquitylome, or acetylome). Use ```panoply_unified_workflow``` to run this pipeline on multiple ome types in parallel. 

This pipeline executes the following modules:

* [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)
* [panoply_normalize_ms_data_report](./Report-Modules%3A-panoply_normalize_ms_data_report)
* [panoply_rna_protein_correlation](./Data-Analysis-Modules%3A-panoply_rna_protein_correlation)
* [panoply_rna_protein_correlation_report](./Report-Modules%3A-panoply_rna_protein_correlation_report)
* [panoply_harmonize](./Support-Modules%3A-panoply_harmonize)
* [panoply_sampleqc](./Data-Preparation-Modules%3A-panoply_sampleqc)
* [panoply_sampleqc_report](./Report-Modules%3A-panoply_sampleqc_report)
* [panoply_cna_setup](./Support-Modules%3A-panoply_cna_setup)
* [panoply_cna_correlation](./Data-Analysis-Modules%3A-panoply_cna_correlation)
* [panoply_cna_correlation_report](./Report-Modules%3A-panoply_cna_correlation_report)
* [panoply_association_workflow](./Data-Analysis-Modules%3A-panoply_association)
* [panoply_accumulate](./Support-Modules%3A-panoply_accumulate)
* [panoply_nmf](./Data-Analysis-Modules%3A-panoply_nmf)
* [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea)
* [panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis) (optional, proteome data only)
* [panoply_download](./Support-Modules%3A-panoply_download)

## Input

Required inputs:

* ```job_identifier```: (String) An identifier name given to the job
* ```ome_type```: (String) Type of omics data being run (ex: "proteome", "phosphoproteome" etc.)
* ```input_pome```: (`.gct` file) Input ome data matrix
* ```input_rna```: (`.gct` file, default = this.rna_ss) Input rna data matrix (optional)
* ```input_cna```: (`.gct` file, default = this.cna_ss) Input cna data matrix (optional)
* ```sample_annotation```: (File, default = this.annotation_ss) The annotation file for the given ome data
* ```groups_file```: (File, default = this.groups_ss) The annotation file for the given ome data, subsetted to annotations of interest and used for analyses. Can be manually overriden for relevant modules (below).
* ```yaml```: (`.yaml` file) parameters in `yaml` format

Optional inputs:
* ```run_ptmsea```: (String, default = ```false```) ```true``` if phosphoproteome data is being run and ptmsea analysis is desired or "false" to skip (this module takes additional time and higher cost to run)
* ```run_cmap```: (String, default = ```false```) ```true``` if proteome data is being run and cmap analysis is desired or ```false``` to skip (this module takes additional time and higher cost to run)
* ```cna_groups```: please see `groupsFile` input parameter in [panoply_cna_setup](./Support-Modules%3A-panoply_cna_setup)
* ```groups_file_association```: please see `groupsFile` input parameter in [panoply_association](./Data-Analysis-Modules%3A-panoply_association). Overrides default `groups_file`.
* ```groups_file_nmf```: please see `groupsFile` input parameter in [panoply_nmf](./Data-Analysis-Modules%3A-panoply_nmf). Overrides default `groups_file`.
* ```normalizeProteomics```: please see `normalizeProteomics` input parameter in [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)

* ```cmap_n_permutations```: please see `n_permutations` input parameter in [panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis)
* ```groups_file_cmap_enrichment```: please see `cmap_enrichment_groups` input parameter in [panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis). Overrides default `groups_file`.
* ```subset_list_file```: please see `subset_list_file` input parameter in [panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis)
* ```cmap_level5_data```: please see `cmap_level5_data` input parameter in [panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis)
* ```annotation_pathway_db```: please see `annotation_pathway_db` input parameter in [panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis)
* ```subset_bucket```: please see `subset_bucket` input parameter in [panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis)

* ```ndigits```: please see `ndigits` input parameter in [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)
* ```na_max```: please see `naMax` input parameter in [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)
* ```sample_na_max```: please see `sample_na_max` input parameter in [panoply_association](./Data-Analysis-Modules%3A-panoply_association)
* ```min_numratio_fraction```: please see `minNumratioFraction` input parameter in [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)
* ```nmiss_factor```: please see `nmiss_factor` input parameter in [panoply_association](./Data-Analysis-Modules%3A-panoply_association)
* ```sd_filter_threshold```: please see `sdFilterThreshold` input parameter in [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)
* ```duplicate_gene_policy```: please see `duplicate_gene_policy` input parameter in [panoply_harmonize](./Support-Modules%3A-panoply_harmonize)
* ```gene_id_col```: please see `gene_id_col` input parameter in [panoply_harmonize](./Support-Modules%3A-panoply_harmonize)
* ```organism```: 

* ```standalone```: (String, default = ```false```) set to ```true``` to run as a self-contained module; when running `panoply_main` pipeline use ```false```
* ```geneset_db```: (String, default = this.gseaDB) please see `gene_set_database` input parameter in [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea)
* ```ptm_db```: (String, default = this.ptmseaDB) this is the `gene_set_database` for running ```ptmsea```


## Output

`panoply_main` produces the follow outputs:

* `panoply_full.tar`:
When opened, contains the following folders:
	- `association`: all results from [panoply_association_workflow](./Data-Analysis-Modules%3A-panoply_association)
	- `so_nmf`: all results from [panoply_nmf](./Data-Analysis-Modules%3A-panoply_nmf)
	- `cna`: all results from [panoply_cna_setup](./Support-Modules%3A-panoply_cna_setup) and [panoply_cna_correlation](./Data-Analysis-Modules%3A-panoply_cna_correlation)
	- `data`: folder containing data needed to run `panoply_main` analysis
	- `harmonized-data`: all results from [panoply_harmonize](./Support-Modules%3A-panoply_harmonize)
	- `normalized-data`: all results from [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)
	- `rna`: all results from [panoply_rna_protein_correlation](./Data-Analysis-Modules%3A-panoply_rna_protein_correlation)
	- `sample-qc`: all results from [panoply_sampleqc](./Data-Preparation-Modules%3A-panoply_sampleqc)
	- `ssgsea_assoc`: all results from running `panoply_ssgsea` after `panoply_association`. Please see [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea) for more information
	- `ssgsea_ome`: all results from running `panoply_ssgsea` on the omics input data. Please see [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea) for more information
	- `ssgsea_rna`: all results from running `panoply_ssgsea` on the rna input data. Please see [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea) for more information
* `summary_and_ssgsea.tar`: contains a smaller version of `panoply_full.tar` with large data files removed for easier download and viewing. The directory structure and folders are the same as listed above but contain fewer data files and/or contain only pdf/png results from the given analysis.
* `cmap_output.tar`(if `panoply_cmap_analysis` is run): contains results for the cmap analysis in the `cmap` folder of the opened `.tar` file. Please see [panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis) for more information on these results.
* `cmap_ssgsea_output.tar`(if `panoply_cmap_analysis` is run): contains results for the cmap annotate ssgsea. Please see [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea) for more information on these results.
* `norm_report.html`: report produced by [panoply_normalize_ms_data_report](./Report-Modules%3A-panoply_normalize_ms_data_report)
* `rna_corr_report.html`: report produced by [panoply_rna_protein_correlation_report](./Report-Modules%3A-panoply_rna_protein_correlation_report)
* `cna_corr_report.html`: report produced by [panoply_cna_correlation_report](./Report-Modules%3A-panoply_cna_correlation_report)
* `sample_qc_report.html`: report produced by [panoply_sampleqc_report](./Report-Modules%3A-panoply_sampleqc_report)
* `association_report.html`: report produced by [panoply_association_report](./Report-Modules%3A-panoply_association_report)
* `so_nmf_report.html`: report produced by [panoply_nmf_report](./Report-Modules%3A-panoply_nmf_report)
* `so_nmf_ssgsea_report.html`: report produced by [panoply_ssgsea_report](./Report-Modules%3A-panoply_ssgsea_report) for NMF clustering results
* `normalized_data_table.gct`: normalized data matrix produced by [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)



