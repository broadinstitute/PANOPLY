# ```panoply_clumps_ptm_postprocess```

## Description

This module wrangles the results of [panoply_clumps_ptm](./Data-Analysis-Modules%3A-panoply_clumps_ptm), creating summary figures and visualizations. It creates a dotplot showing the top proteins with clusters of differentially-upregulated, and differentially down-regulated PTM-sites, for each annotation value tested. It also creates PyMol figures for the top significant proteins, with the location of modified residues colored.

## Input

### Required inputs:

* ```results_tar```: (`.tar` file) tar file containing the results of [panoply_clumps_ptm](./Data-Analysis-Modules%3A-panoply_clumps_ptm)
* ```output_prefix```: (String, default="results") prefix used to name the output files
* ```yaml_file```: (`.yaml` file) master-parameters.yaml

### Optional inputs:

* ```fdr_threshold```: (Float) FDR threshold for significance
* ```pymol_gen```: (Boolean, default=`true`) if `TRUE`, PyMol figures will be generated
* ```pymol_upper_limit```: (Int, default=10) maximum number of PyMol figures that should be generated per test


## Output

* ```results```: (`"${output_prefix}_clumps_ptm_full_results.tar"`) all results and figures of ClumpsPTM, including data-tables and raw results
* ```figures_only```: (`"${output_prefix}_clumps_figures.tar"`) summary figures and PyMol only 