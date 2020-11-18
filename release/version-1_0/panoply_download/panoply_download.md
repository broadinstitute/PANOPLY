Documentation at https://github.com/broadinstitute/PANOPLY/blob/release-1_0/release/version-1_0/panoply_download/panoply_download.md

# ```panoply_download```

## Description
Assembles the outputs from various analysis and report modules in the `panoply_main` pipeline into a coherent directory structure and creates two final output tarballs containing: (i) reports and relevant output tables/plots for all analysis modules and (ii) the complete results from the entire pipeline.

## Input

Required inputs:
* ```cons_clust_tar```: (`.tar` file) tarball from panoply_cons_clust
* ```ssgsea_ome_tar```: (`.tar` file)
* ```ssgsea_rna_tar```: (`.tar` file)
* ```analysisDir```: (String) name of analysis directory
* ```ssgsea_assoc_tars```: (Array[File]) array of `.tar` files from ssGSEA results from panoply_association
* ```ssgsea_clust_tars```: (Array[File]) array of `.tar` files from ssGSEA results from panoply_cons_clust
* ```output_prefix```: (String) prefix added to the summary_tar and full_tar output names (usually set to omics type eg: 'proteome' etc.)

Optional inputs:
* ```ptmsea```: (File)
* ```summary_tar```: (String, default = 'panoply_main_summary.tar') name of output summary `.tar` file
* ```full_tar```: (String, default = 'panoply_main_full.tar') name of output full `.tar` file
* ```ssgsea_assoc_dir```: (String, default = 'ssgsea_assoc') directory containing ssGSEA Association results
* ```ssgsea_clust_dir```: (String, default = 'ssgsea_clust') directory containing ssGSEA clustering results

## Output

* Summary `.tar` file
* Full `.tar` file