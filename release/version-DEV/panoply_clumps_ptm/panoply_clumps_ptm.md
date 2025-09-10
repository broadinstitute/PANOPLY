# ```panoply_clumps_ptm```

## Description

This module runs [Clumps-PTM](https://github.com/getzlab/CLUMPS-PTM/tree/main), a top-down spatial-proteomics analysis tool that identifies proteins with nearby-clusters of differentially-regulated PTM-sites (phosphorylation, acetylation, and/or ubiquitination). The algorithm was adapted from the CLUMPS method for detecting clusters of mutations in 3D protein structures; it calculates a weighted average proximity score across all differentially-modified residue pairs in a given protein, with weights given according to logFC and significance. An empirical p-value is calculated by permuting across the possible PTM-sites within the protein, before correction for multiple-testing. A full description of the algorithm can be found in the Method Details of [Geffen et al. 2023](https://doi.org/10.1016/j.cell.2023.07.013).

## Input

### Required inputs:

* ```diff_exp_file```: (`.tsv` file) results file from [panoply_clumps_ptm_diffexp](./Data-Analysis-Modules%3A-panoply_clumps_ptm_diffexp), containing differential expression results for all PTM -omes, for a given annotation
* ```var_sites_file```: (`.tsv` file) filtered mapping file (`filt_results`) from [panoply_clumps_ptm_mapping](./Data-Preparation-Modules%3A-panoply_clumps_ptm_mapping), containing all varaible sites with valid PDB coordinates

* ```PDB_ref_bucket```: (String) Google-Cloud Bucket containing a tarred copy of the PDB structural archive (i.e. `https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/`). A public bucket, pulled from a frozen 2025 snapshot, can be found at: `"gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/pdbs_2025_tars/"`
	* ```PDB_DIR```: Internal parameter listing the files to import from `PDB_ref_bucket`

* ```output_prefix```: (String, default="results") prefix used to name the output tar file
* ```yaml_file```: (`.yaml` file) master-parameters.yaml

## Optional inputs:

* ```run_combined```: (Boolean, default=`true`) if `TRUE` analysis will be run on all PTM datasets combined, in addition to each -ome separately
* ```weight_col```: (String, default="logFC") column from differential-expression dataset to use as weights in ClumpsPTM

* ```accession_col```: (String, default="description") GCT rdesc column with protein accession IDs; must use the same ID type as the provided `FASTA_ref_file` file.
* ```variable_sites_col```: (String, default="variableSites") GCT rdesc column with PTM variable site(s) (e.g. 'T527t')

* ```DEBUG_MODE```: (Boolean, default=`false`) Debugging toggle; if `true`, a small subset of proteins will be analyzed. **Should be turned off for analysis.**

## Output

* ```results```: (`.tar` file) 



## References

1. Geffen, Y. et al. **Pan-cancer analysis of post-translational modifications reveals shared patterns of protein regulation.** Cell 186, 3945-3967.e26 (2023).
