# ```panoply_preprocess_gct```

## Description

This module collapses a feature-level GCT to a gene-centric (for ssGSEA) or site-centric (for PTM-SEA) level, as an appropriate input GCT for the [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea) module. 

## Input

### Required inputs:

* ```input_ds```: (`.gct` file) input GCT file
* ```yaml_file```: (`.yaml` file) master-parameters.yaml
* ```output_prefix```: (String) File prefix for output files.

### Optional inputs:
* ```acc_type```: (String) Type of accession number in 'rid' object in GCT file ("uniprot", "refseq" (default), "symbol").
* ```id_type```: (String) Notation of site-ids: 'sm' - Spectrum Mill (default); 'wg' - Web Gestalt; 'ph' - Philosopher. Only relevant for PTM-SEA.
* ```id_type_out``` (String) Type of site id for output: 'uniprot'(default), 'refseq', 'seqwin'. Only relevant for PTM-SEA.
* ```level``` (String) Mode of report:
  + 'ssc' - single-site-centric
  + 'gc' - gene-centric (default)
  + 'gcr' - gene-centric-redundant
* ```loc``` (Boolean) If TRUE only fully localized sites will be considered (default: TRUE). Localization infromation is expected to be encoded in the site identifier. Respective parsing rules are determined by '--id_type'.
* ```gene_col```: (String) Name of column listing gene names; used for gene centric reports (default: "geneSymbol").
* ```humanize_gene```: (Boolean) If TRUE, gene symbols will be capitalized; can be used to crudely humanize mouse or rat gene symbols.
* ```seqwin_col```: (String) "Column containing flanking sequences, separated by '|'. Only relevant for PTM-SEA and if '--id_type_out' = 'seqwin' (default: 'VMsiteFlanks').
* ```SGT_col```: (String) Column used to collpase subgroup-top (SGT) reports (default: "subgroupNum). Only relevant for Spectrum Mill protein reports.
* ```mod_res```: (String) Modified residues, e.g. "S|T|Y" or "K" (default: "S|T|Y").
* ```mod_type```: (String) Type of post-translational modification, e.g "p" for phospho (default) or "ac" for acetylation
* ```mode```: (String) Determines how multiple features (e.g. proteins, PTM sites, etc.) mapping to the same gene symbol will be aggregated: 
  + "mean" - mean
  + "median" - median
  + "sd - most variable (standard deviation) across sample columns
  + "SGT" - subgroup top: first subgroup in protein group (Spectrum Mill)
  + "abs.max" - for log-transformed, signed p-values"

## Output

* ```result```: (`.gct` file) Preprocessed GCT file, appropriate for use in ssGSEA or PTM-SEA
