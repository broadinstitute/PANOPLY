# ```panoply_unified_assemble_results```

## Description
Performs final assembly of all pipeline analysis results and reports from `panoply_unified_workflow`. 

## Input

Required inputs: 
* ```output_results_zip```: (String, default = "all_results.zip") results output `.zip` file name 
* ```output_reports_zip```: (String, default = "all_reports.zip") reports output `.zip` file name

Optional inputs:
* ```main_full```: (Array[File]) array of full `.tar` results files from the panoply_main proteogenomic analysis 
* ```main_summary```: (Array[File]) array of summary `.tar` results files from the panoply_main proteogenomic analysis
* ```norm_report```: (Array[File]) array of `.html` report files from `panoply_normalize_ms_data_report` in the `panoply_main` proteogenomic analysis
* ```rna_corr_report```: (Array[File]) array of `.html` report files from `panoply_rna_protein_correlation_report` in the `panoply_main` proteogenomic analysis
* ```cna_corr_report```: (Array[File]) array of `.html` report files from `panoply_cna_correlation_report` in the `panoply_main` proteogenomic analysis
* ```sampleqc_report```: (Array[File]) array of `.html` report files from `panoply_sampleqc_report` in the `panoply_main` proteogenomic analysis
* ```assoc_report```: (Array[File]) array of `.html` report files from `panoply_association_report` in the `panoply_main` proteogenomic analysis
* ```cons_clust_report```: (Array[File]) array of `.html` report files from `panoply_cons_clust_report` in the `panoply_main` proteogenomic analysis
* ```blacksheep_tar```: (Array[File]) array of `.tar` results files from `panoply_blacksheep` outlier analysis
* ```blacksheep_report```: (Array[File]) array of `.html` report files from `panoply_blacksheep_report` 
* ```cmap_output```: (Array[File]) array of `.tar` output files from `panoply_cmap_analysis` in the `panoply_main` proteogenomic analysis
* ```cmap_ssgsea_output```: (Array[File]) array of `.tar` output files from running ssGSEA on `panoply_cmap_analysis` in the `panoply_main` proteogenomic analysis
* ```nmf_results```: (File) `.tar` file with all results from `panoply_nmf_workflow`
* ```nmf_report```: (File) `.tar` file with all reports from `panoply_nmf_workflow`
* ```immune_tar```: (File) `.tar` file results from `panoply_immune_analysis`
* ```immune_report```: (File) `.html` file showing results from `panoply_immune_analysis` analysis

## Output

* all_results `.zip` file
* all_reports `.zip` file