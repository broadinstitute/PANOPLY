# ```panoply_main_internal```

## Description

Performs single-omic analysis for a given -ome (e.g. proteome, transcriptomic, phosphoproteome, etc.). Use [panoply_main](./Pipelines%3A-panoply_main) to run this pipeline along with additional proteogenomic analyses, if RNA and CNA data are available.

This pipeline executes the following modules:

* [panoply_association_workflow](./Workflows%3A-panoply_association_workflow)
* [panoply_blacksheep_workflow](./Workflows%3A-panoply_blacksheep_workflow)
* [panoply_ssgsea_workflow](./Workflows%3A-panoply_ssgsea_workflow)
* [panoply_nmf_internal_workflow](./Workflows%3A-panoply_nmf_internal_workflow)
* [panoply_immune_analysis_workflow](./Workflows%3A-panoply_immune_analysis_workflow) (transcriptomics only)
* [single-sample PTM-SEA](./Workflows%3A-panoply_ssgsea_workflow) (phosphoproteome only)

## Input

### Required inputs:

* ```job_identifier```: (String) An identifier name given to the job
* ```input_ome```: (`.gct` file) Input ome data matrix
* ```sample_annotation```: (File, default = this.annotation_ss) The annotation file for the given ome data
* ```groups_file```: (File, default = this.groups_ss) The annotation file for the given ome data, subsetted to annotations of interest and used for analyses. Can be manually overriden for relevant modules (below).
* ```yaml```: (`.yaml` file) parameters in `yaml` format

### Module Toggles:
* ```run_ptmsea```: (String, default = ```false```) ```true``` if phosphoproteome data is being run and ptmsea analysis is desired or "false" to skip (this module takes additional time and higher cost to run)
* ```run_nmf```: (String, default = ```true```) ```true``` if nmf analysis is desired or ```false``` to skip

### Optional inputs:

* ```sample_na_max```: please see `sample_na_max` input parameter in [panoply_association](./Data-Analysis-Modules%3A-panoply_association)
* ```nmiss_factor```: please see `nmiss_factor` input parameter in [panoply_association](./Data-Analysis-Modules%3A-panoply_association)
* ```duplicate_gene_policy```: please see `duplicate_gene_policy` input parameter in [panoply_harmonize](./Support-Modules%3A-panoply_harmonize)
* ```gene_id_col```: please see `gene_id_col` input parameter in [panoply_harmonize](./Support-Modules%3A-panoply_harmonize)

* ```standalone```: (String, default = ```false```) set to ```true``` to run as a self-contained module; when running `panoply_main` pipeline use ```false```
* ```geneset_db```: (String, default = this.gseaDB) please see `gene_set_database` input parameter in [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea)
* ```ptm_db```: (String, default = this.ptmseaDB) this is the `gene_set_database` for running ```ptmsea```


#### Groups File Overrides
By default, all analyses take `groups_file` as the default annotations-file. However, for certain modules it may be desireable to exclude, add, or otherwise customize the annotations used. The following parameters allow groups-files to be overridden on a module-specific basis:

* ```groups_file_association```: please see `groupsFile` input parameter in [panoply_association](./Data-Analysis-Modules%3A-panoply_association). Overrides default `groups_file`.
* ```groups_file_blacksheep```: please see `groupsFile` input parameter in [panoply_blacksheep](./Data-Analysis-Modules%3A-panoply_blacksheep). Overrides default `groups_file`.
* ```groups_file_immune```: please see `groupsFile` input parameter in [panoply_immune_analysis](./Data-Analysis-Modules%3A-panoply_immune_analysis). Overrides default `groups_file`.
* ```groups_file_nmf```: please see `groups_file` input parameter in [panoply_nmf_postprocess](./Data-Analysis-Modules%3A-panoply_nmf_postprocess). Overrides default `groups_file`.



## Output

`panoply_main_internal` produces a tarfile and report file for each of its modules

* ```association_tar```: tarfile output of [panoply_association](./Data-Analysis-Modules%3A-panoply_association)
* ```ssgsea_assoc_tars```: ssGSEA results for all annotations analyzed in [panoply_association](./Data-Analysis-Modules%3A-panoply_association)
* ```association_report```: report produced by [panoply_association_report](./Report-Modules%3A-panoply_association_report)

* ```blacksheep_tar```: tarfile output of [panoply_blacksheep](./Data-Analysis-Modules%3A-panoply_blacksheep)
* ```blacksheep_report```: report produced by [panoply_blacksheep_report](./Report-Modules%3A-panoply_blacksheep_report)

* ```ssgsea_ome_tar```: tarfile output of [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea)
* ```ssgsea_ome_report```: report produced by [panoply_ssgsea_report](./Report-Modules%3A-panoply_ssgsea_report)

* ```so_nmf_results```: tarfile output of [panoply_nmf](./Data-Analysis-Modules%3A-panoply_nmf), containing the results of NMF clustering
* ```so_nmf_figures```: tarfile output of [panoply_nmf_postprocess](./Data-Analysis-Modules%3A-panoply_nmf_postprocess), containing various postprocessing and figures from NMF analysis
* ```so_nmf_report```: report produced by [panoply_nmf_report](./Report-Modules%3A-panoply_nmf_report)
* ```so_nmf_ssgsea_tar```: ssGSEA results on the W-matrix from NMF results; GCT produced by [panoply_nmf_postprocess](./Data-Analysis-Modules%3A-panoply_nmf_postprocess)
* ```so_nmf_ssgsea_report```: report produced by [panoply_ssgsea_report](./Report-Modules%3A-panoply_ssgsea_report) for NMF clustering results

* `normalized_data_table.gct`: normalized data matrix produced by [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)

* ```immune_analysis_tar```: tarfile output of [panoply_immune_analysis](./Data-Analysis-Modules%3A-panoply_immune_analysis), if applicable
* ```immune_report```: report produced by [panoply_immune_analysis_report](./Report-Modules%3A-panoply_immune_analysis_report), if applicable

* ```ptmsea_tar```: tarfile output of PTM-SEA results produced by [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea), if applicable
* ```ptmsea_ome_report```: report of PTM-SEA results produced by [panoply_ssgsea_report](./Report-Modules%3A-panoply_ssgsea_report), if applicable




