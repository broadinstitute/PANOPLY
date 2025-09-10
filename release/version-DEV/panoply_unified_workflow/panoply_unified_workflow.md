# ```panoply_unified_workflow```

## Description
Performs proteogenomic analysis on multiple omics types in parallel (ie ```panoply_main```) and leverages additional analysis modules such as ```panoply_nmf```, ```panoply_clumps_ptm_workflow```, and ```panoply_metaboanalyst```. It can optionally be used to normalize and filter -omic datasets using ```panoply_normalize_filter_workflow```.

This pipeline executes the following modules and pipelines:

* [panoply_normalize_filter_workflow](./Workflows%3A-panoply_normalize_filter_workflow)
* [panoply_main](./Pipelines%3A-panoply_main)
* [panoply_clumps_ptm_workflow](./Workflows%3A-panoply_clumps_ptm_workflow)
* [panoply_nmf_workflow](./Workflows%3A-panoply_nmf_workflow)
* [panoply_metaboanalyst](./Workflows%3A-panoply_metaboanalyst_workflow), if applicable

## Input

### Required inputs:

* ```yaml```: (`.yaml` file) parameters in `yaml` format
* ```job_id```: (String) An identifier name given to the job
* ```groups_file```: (File, default = this.groups_ss) The annotation file for the given ome data, subsetted to annotations of interest and used for analyses. Can be manually overriden for relevant modules (below).

#### -Omic Inputs
**NOTE:** at least one *must* be specified):

* ```prote_ome```: (`.gct` file, default = this.proteome_ss) proteome data matrix
* ```phospho_ome```: (`.gct` file, default = this.phosphoproteome_ss) phosphoproteome data matrix
* ```acetyl_ome```: (`.gct` file, default = this.acetylome_ss) acetylome data matrix
* ```ubiquityl_ome```: (`.gct` file, default = this.ubiquitylome_ss) ubiquitylome data matrix
* ```nglyco_ome```: (`.gct` file, default = this.nglycoproteome_ss) nglycoproteome data matrix
* ```methyl_ome```: (`.gct` file, default = this.methylation_ss) methylation data matrix

#### Non-Proteomic Inputs
* Genomic Inputs: must specify **both** RNA and CNA data, in order to run proteogenomic analyses in [panoply_main](./Pipelines%3A-panoply_main) 
	* ```rna_data```: (`.gct` file, default = this.rna_v3_ss) Input rna data matrix
	* ```cna_data```: (`.gct` file, default = this.cna_ss) Input cna data matrix
* ```metabol_ome```: (`.gct` file, default = this.metabolome_ss) metabolome data matrix

#### Module Toggles:
* ```run_ptmsea```: (String, default = ```false```) ```true``` if phosphoproteome data is being run and ptmsea analysis is desired or "false" to skip (this module takes additional time and higher cost to run)
* ```run_clumps```: (String, default = ```false```) ```true``` if PTM data is being run and ClumpsPTM analysis is desired, or "false" to skip (this module takes additional time and higher cost to run)
* ```run_cmap```: (String, default = ```false```) ```true``` if proteome data is being run and cmap analysis is desired or ```false``` to skip (this module takes additional time and higher cost to run)
* ```run_nmf```: (String, default = ```true```) ```true``` if nmf analysis is desired or ```false``` to skip
* ```run_omicsev```: (String, default = ```true```) ```true``` if OmicsEV analysis is desired or ```false``` to skip


### Optional Inputs

* ```normalizeProteomics```: please see `normalizeProteomics` input parameter in [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)
* ```filterProteomics```: please see `filterProteomics` input parameter in [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_filter)

#### Groups File Overrides
By default, all analyses take `groups_file` as the default annotations-file. However, for certain modules it may be desireable to exclude, add, or otherwise customize the annotations used. The following parameters allow groups-files to be overridden on a module-specific basis:

* ```groups_file_clumpsptm```: please see `groupsFile` input parameter in [panoply_clumps_ptm_diffexp](./Data-Analysis-Modules%3A-panoply_clumps_ptm_diffexp). Overrides default `groups_file`.
* ```groups_file_metaboanlayst```: please see `groups_file` input parameter in [panoply_metaboanalyst](./Data-Analysis-Modules%3A-panoply_metaboanalyst). Overrides default `groups_file`.
* ```groups_file_nmf```: please see `groups_file` input parameter in [panoply_nmf_postprocess](./Data-Analysis-Modules%3A-panoply_nmf_postprocess). Overrides default `groups_file`.


## Output
This pipeline outputs two `.zip	` files:

1. `all_results.zip` 
This file contains complete results from all pipelines and modules run. The directory structure and results are formatted as follows:

* `results`:
	- `proteogenomic_analysis`: contains all results from `panoply_main` including a folder called `all_html_reports` which contains all reports produced from the `panoply_main` pipeline. If `panoply_cmap_analysis` was run, a folder named `proteome_cmap_analysis` will be present containing results from this module. Please see the outputs from [panoply_main](./Pipelines%3A-panoply_main) for more information.
	- `rna`: contains all results from [panoply_main_internal](./Pipelines%3A-panoply_main_internal) run on RNA data. Please see the outputs from [panoply_main_internal](./Pipelines%3A-panoply_main_internal) for more information.
	- `clumps_ptm`: contains results from [panoply_clumps_ptm_workflow](./Workflows%3A-panoply_clumps_ptm_workflow), if applicable
	- `metaboanalyst`: contains results from [panoply_metaboanalyst](./Workflows%3A-panoply_metaboanalyst_workflow), if applicable
	- `nmf`: contains results from [panoply_nmf](./Workflows%3A-panoply_nmf_workflow)

2. `all_reports.zip`
This file contains all reports produced from every pipeline and module run. The directory structure for this file is formatted as follows:

* `reports`:
	- `proteogenomic_analysis`: contains all reports generated by the modules in [panoply_main](./Pipelines%3A-panoply_main)
	- `rna`: contains all reports generated by the modules in [panoply_main_internal](./Pipelines%3A-panoply_main_internal) on RNA data
	- `clumps_ptm`: contains all reports generated by [panoply_clumps_ptm_report](./Report-Modules%3A-panoply_clumps_ptm_report)
	- `metaboanalyst`: contains all reports generated by [panoply_metaboanalyst_report](./Report-Modules%3A-panoply_metaboanalyst_report)
	- `nmf`: contains all reports generated by [panoply_nmf_report](./Report-Modules%3A-panoply_nmf_report) and [panoply_ssgsea_report](./Report-Modules%3A-panoply_ssgsea_report)
