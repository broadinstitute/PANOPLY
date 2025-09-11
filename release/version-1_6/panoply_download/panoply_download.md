# ```panoply_download```

## Description
Assembles the outputs from various analysis and report modules in the `panoply_main` pipeline into a coherent directory structure and creates two final output tarballs containing: (i) reports and relevant output tables/plots for all analysis modules and (ii) the complete results from the entire pipeline.

## Input

### Required inputs:

* ```analysisDir```: (String) name of analysis directory
* ```output_prefix```: (String) prefix added to the summary_tar and full_tar output names (usually set to omics type eg: 'proteome' etc.)

* ```association_tar```: (`.tar` file) results from [panoply_association](./Data-Analysis-Modules%3A-panoply_association)
* ```blacksheep_tar```: (`.tar` file) tarball from [blacksheep outlier analysis](./Data-Analysis-Modules%3A-panoply_blacksheep)
* ```ssgsea_assoc_tars```: (Array[File]) array of `.tar` files from ssGSEA results from [panoply_association](./Workflows%3A-panoply_association_workflow)
* ```ssgsea_ome_tar```: (`.tar` file) single-sample GSEA results of [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea)


### Optional inputs:

* ```summary_tar```: (String, default = 'panoply_main_summary.tar') name of output summary `.tar` file
* ```full_tar```: (String, default = 'panoply_main_full.tar') name of output full `.tar` file
* ```ssgsea_assoc_dir```: (String, default = 'ssgsea_assoc') directory containing ssGSEA Association results

#### Optional Analysis Modules

* ```cna_corr_tar```: (`.tar` file) tarball containing all proteogenomic analysis results from [panoply_main](./Pipelines%3A-panoply_main)
* ```so_nmf_results```: (`.tar` file) tarball from single-omic [nmf clustering](./Data-Analysis-Modules%3A-panoply_nmf)
* ```so_nmf_figures```: (`.tar` file) tarball with figures from [nmf postprocessing](./Data-Analysis-Modules%3A-panoply_nmf_postprocess)
* ```so_nmf_ssgsea_tar```: (`.tar` file) tarball with ssgsea results on the W-Matrix from single-omic NMF clustering

#### Optional QC Modules
* ```omicsev_tar```: (`.tar` file) results from [panoply_omicsev](./Data-Analysis-Modules%3A-panoply_omicsev)
* ```cosmo_tar```: (`.tar` file) results from [panoply_cosmo](./Data-Preparation-Modules%3A-panoply_cosmo)



#### -ome Specific Modules
* ```immune_analysis_tar```: (File) RNA-only; immune analysis results from [panoply_immune_analysis_workflow](./Workflows%3A-panoply_immune_analysis_workflow)
* ```ptmsea```: (File) phosphoproteome-only; PTM-SEA results from [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea)

## Output

* ```summary```: Summary `.tar` file
* ```full```: Full `.tar` file