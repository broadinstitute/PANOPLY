
# ```panoply_nmf_assemble_results```

## Description

This module compiles the outputs of multiple [panoply_nmf_internal_workflow](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) workflow-runs into a single results tar and a single reports tar.

## Input

* Multi-omic NMF results (if run)
	- ```mo_nmf_results```: (File) Tar file containing results of [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf)
	- ```mo_nmf_figures```: (File) Tar file containing figures and analyses from [panoply_nmf_postprocess](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf_postprocess)
	- ```mo_nmf_report```: (File) Summary report of NMF clustering pipeline.
	- ```mo_nmf_ssgsea_tar```: (File) Results of [panoply_ssgsea](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea) applied to the NMF results.
	- ```mo_nmf_ssgsea_report```: (File) Summary report of ssGSEA applied to the NMF results.


* Single-omic NMF results (if run)
	- ```so_nmf_results```: (Array[File]) Tar files containing results of [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf)
	- ```so_nmf_figures```: (Array[File]) Tar files containing figures and analyses from [panoply_nmf_postprocess](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf_postprocess)
	- ```so_nmf_report```: (Array[File]) Summary reports of NMF clustering pipeline.
	- ```so_nmf_ssgsea_tar```: (Array[File]) Results of [panoply_ssgsea](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea) applied to the NMF results.
	- ```so_nmf_ssgsea_report```: (Array[File]) Summary reports of ssGSEA applied to the NMF results.

* Sankey Results (if run)
	- ```sankey_tar```: (File) Tar file containing sankey diagram results of [panoply_sankey](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_sankey)
	- ```sankey_report```: (File) Summary reports of sankey-diagram comparisons


## Output

* ```nmf_resuls.tar```: Tarred results of NMF clustering and GSEA analysis, from all NMF runs (single-omic and multi-omic). Also contains sankey diagrams (if run).
* ```nmf_reports.tar```: Tarred reports for NMF clustering and GSEA analysis, from all NMF runs (single-omic and multi-omic). Also contains sankey diagrams (if run).