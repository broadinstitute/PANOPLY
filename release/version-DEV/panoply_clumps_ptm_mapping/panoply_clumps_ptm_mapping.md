# ```panoply_clumps_ptm_mapping```

## Description

This module takes input PTM datasets and maps each feature's variable site to atomic coordinates.

While mapping using AlphaFold structures is possible using the [ClumpsPTM](https://github.com/getzlab/CLUMPS-PTM/tree/main) toolset, but is not currently integrated into this module.

## Input

### Required inputs:

* ```pSTY_gct```: (`.gct` file) phosphoproteome data matrix
* ```acK_gct```: (`.gct` file) acetylome data matrix
* ```ubK_gct```: (`.gct` file) ubiquitylome data matrix

* ```output_prefix```: (String, default="results") prefix used to name the output tar file
* ```yaml_file```: (`.yaml` file) master-parameters.yaml


#### Mapping Databases 

* ```PDB_ref_bucket```: (String) Google-Cloud Bucket containing a tarred copy of the PDB structural archive (i.e. `https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/`). A public bucket, pulled from a frozen 2025 snapshot, can be found at: `"gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/pdbs_2025_tars/"`
	* ```PDB_DIR```: Internal parameter listing the files to import from `PDB_ref_bucket`
* ```FASTA_ref_file```: (File) Reference FASTA file with all relevant sequences for your dataset. Used to BLAST against Uniprot sequences. 
* ```UNIPROT_SWISSPROT```: (File) Reference FASTA file with all relevant UNIPROT sequences, to which your sequences will be BLASTed.
* ```SIFTS_DB```: (File) SIFTS database containing mapping between UNIPROT IDs and PDB IDs.

### Optional inputs:

* ```FASTA_sep_type```: (String) separator for FASTA sequences in FASTA reference file. Supported values are "cptac" for ' ', or "gencode" for "|"
* ```accession_col```: (String, default="description") GCT rdesc column with protein accession IDs; must use the same ID type as the provided `FASTA_ref_file` file.
* ```gene_column```: (String, default="geneSymbol") GCT rdesc column with HUGO gene symbols
* ```variable_sites_col```: (String, default="variableSites") GCT rdesc column with PTM variable site(s) (e.g. 'T527t')
* ```variable_sites_sep```: (String, default=" ") separator for variable sites (e.g. ' ' is the separator for 'T972t S977s')
* ```keep_multi_sites```: (Boolean, default=`true`) should multi-site PTMs be mapped and analyzed?
* ```filter_duplicate_sites```: (Boolean, default=`true`) should multi-site PTMs that were also observed as single-sites be filtered out?

* ```DEBUG_MODE```: (Boolean, default=`false`) Debugging toggle; if `true`, a small subset of FASTA files will be processed. **Should be turned off for analysis.**

## Output

* ```full_results```: (`.tsv` file, `"${output_prefix}_full_mapped_sites_to_pdbs.tsv"`) TSV file containing all variable sites, including those without valid PDB mappings
* ```filt_results```: (`.tsv` file, `"${output_prefix}_mapped_sites_to_pdbs.tsv"`) TSV file containing only variable sites with valid PDB mappings
* ```var_sites_file```: (`.tsv` file, `"${output_prefix}_var_sites_combined.tsv"`) concatenated table of all variable sites from every input-PTM database, without mappings
* ```mapping_params```: (`.yaml` file) parameters file containing parameters used for mapping; used in [panoply_clumps_ptm_report](./Report-Modules%3A-panoply_clumps_ptm_report.md) to display mapping parameters in the report.


