# ```panoply_clumps_ptm_workflow```

## Description

Workflow adaptation of the [Clumps-PTM](https://github.com/getzlab/CLUMPS-PTM/tree/main) tool, which is a top-down spatial-proteomics analysis tool that identifies proteins with nearby-clusters of differentially regulated PTM-sites. The workflow first generates [differential-expression results](./Data-Analysis-Modules%3A-panoply_clumps_ptm_diffexp) for every feature and [maps PTM variable-sites](./Data-Preparation-Modules%3A-panoply_clumps_ptm_mapping) to atomic coordinates, after which it performs [Clumps-PTM](./Data-Analysis-Modules%3A-panoply_clumps_ptm) analysis to identify proteins with clusters of differentially-regulated sites. After analysis is completed the [panoply_clumps_ptm_postprocess](./Data-Analysis-Modules%3A-panoply_clumps_ptm_postprocess) generates summary figures and PyMol protein-structures, and [panoply_clumps_ptm_report](./Report-Modules%3A-panoply_clumps_ptm_report.md) generates an interactive report.

This workflow is time-intensive and can become prohibitively expensive if too many annotations are chosen; it is recommended to select relatively few annotations when providing a `groupsFile`. Additionally, once [panoply_mapping](./Data-Preparation-Modules%3A-panoply_clumps_ptm_mapping) has been run once for a given dataset, it can be skipped in future runs by directly providing the output file in `mapping_file` (assuming the features' accession-numbers and variable-sites IDs have not changed).

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_clumps_ptm_diffexp</code>](./Data-Analysis-Modules%3A-panoply_clumps_ptm_diffexp)         |  performs differential-expression analysis on provided PTM datasets for chosen annotations  |
| [<code>panoply_clumps_ptm_mapping</code>](./Data-Preparation-Modules%3A-panoply_clumps_ptm_mapping)         |  maps PTM sites to atomic coordinates  |
| [<code>panoply_clumps_ptm</code>](./Data-Analysis-Modules%3A-panoply_clumps_ptm)         |  runs ClumpsPTM algorithm  |
| [<code>panoply_clumps_ptm_postprocess</code>](./Data-Analysis-Modules%3A-panoply_clumps_ptm_postprocess)         |  generates figures for ClumpsPTM analysis results  |
| [<code>panoply_clumps_ptm_report</code>](./Report-Modules%3A-panoply_clumps_ptm_report.md)   |  creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the clumps_ptm results |


## Input

### Required inputs:

* ```pSTY_gct```: (`.gct` file) phosphoproteome data matrix
* ```acK_gct```: (`.gct` file) acetylome data matrix
* ```ubK_gct```: (`.gct` file) ubiquitylome data matrix

* ```groupsFile```: (`.csv` file) annotation file, subsetted to annotations of interest for this analysis

* ```output_prefix```: (String, default="results") prefix used to name the output tar file
* ```yaml_file```: (`.yaml` file) master-parameters.yaml


#### Mapping Databases 

* ```PDB_ref_bucket```: (String, default=`"gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/pdbs_2025_tars/"`) Google-Cloud Bucket containing a tarred copy of the PDB structural archive (i.e. `https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/`)
* ```UNIPROT_SWISSPROT```: (File, default=`"gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/reference_files/uniprot_sprot.fasta"`) Reference FASTA file with all relevant UNIPROT sequences, to which your sequences will be BLASTed.
* ```SIFTS_DB```: (File, default=`"gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/reference_files/pdb_chain_uniprot.tsv"`) SIFTS database containing mapping between UNIPROT IDs and PDB IDs.


### Optional inputs:

* ```accession_col```: (String, default="description") GCT rdesc column with protein accession IDs; must use the same ID type as the provided `FASTA_ref_file` file.
* ```variable_sites_col```: (String, default="variableSites") GCT rdesc column with PTM variable site(s) (e.g. 'T527t')

* ```mapping_file```: (`.tsv` file) output file with mapping results from [panoply_mapping](./Data-Preparation-Modules%3A-panoply_clumps_ptm_mapping). If provided directly, `panoply_mapping` will be skipped; useful when re-analyzing a dataset that has previously been mapped.
* ```mapping_params```: (`.yaml` file) parameters file from [panoply_clumps_ptm_mapping](./Data-Preparation-Modules%3A-panoply_clumps_ptm_mapping), containing parameters used for PTM-mapping; allows mapping parameters to be listed in the report.


## Output

`panoply_clumps_ptm_workflow` produces the follow outputs:

* ```clumps_ptm_tar```: (`.tar`) An output tar file that contains the clumps_ptm analysis results 
* ```clumps_ptm_report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
