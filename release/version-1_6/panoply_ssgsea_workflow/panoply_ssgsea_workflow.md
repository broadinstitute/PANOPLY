# ```panoply_ssgsea_workflow```

## Description

This workflow performs single sample Gene Set Enrichment Analysis (ssGSEA) or PTM-Signature Enrichment Analysis (PTM-SEA) on each column of the input data matrix ([panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea)). The module also generates an [interactive report](./Report-Modules%3A-panoply_ssgsea_report.md). Optionally, data can be preprocessed with [panoply_preprocess_gct](./Data-Preparation-Modules%3A-panoply_preprocess_gct) to collapse from protein/peptide-level to gene-centric or site-centric. This workflow executes the following modules:

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_preprocess_gct</code>](./Data-Preparation-Modules%3A-panoply_preprocess_gct)      |  creates a gene-centric or site-centric GCT file suitable for ssGSEA or PTM-SEA  |
| [<code>panoply_ssgsea</code>](./Data-Analysis-Modules%3A-panoply_ssgsea)         |  performs single sample Gene Set Enrichment Analysis (ssGSEA) or PTM-Signature Enrichment Analysis (PTM-SEA) |
| [<code>panoply_ssgsea_report</code>](./Report-Modules%3A-panoply_ssgsea_report.md)   |  creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the ssgsea results |


### panoply_preprocess_gct

This module collapses a feature-level GCT to a gene-centric (for ssGSEA) or site-centric (for PTM-SEA) level GCT. It can optionally be run before [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea) module by toggling the `preprocess_gct` parameter, to collapse a GCT with redundant features into an appropriately-formatted input file.

### panoply_ssgsea

This module performs single sample Gene Set Enrichment Analysis (ssGSEA) or PTM-Signature Enrichment Analysis (PTM-SEA) [1] on each column of the input data matrix. This module is based on the implementation available at the [ssGSEA2.0 GitHub repository](https://github.com/broadinstitute/ssGSEA2.0). 

This is an updated version of the original ssGSEA [2,3] R-implementation. Depending on the input dataset and chosen database (gene sets or PTM signatures), the software performs either ssGSEA or PTM-SEA, respectively. The Molecular Signatures Database ([MSigDB](http://software.broadinstitute.org/gsea/msigdb/)) [4] provides a large collection of curated gene sets.  Gene sets are stored as plain text in  [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) format. A current version of MSigDB gene set collections can be found in the ```db/msigdb``` subfolder. MSigDB gene sets are realeased under [Creative Commons Attribution 4.0 International License](http://software.broadinstitute.org/gsea/msigdb_license_terms.jsp). The license terms can be found in the```db/msigdb``` folder.

File formats supported by ssGSEA2.0/PTM-SEA are Gene Cluster Text [GCT v1.2](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29) or [GCT v1.3](https://clue.io/connectopedia/gct_format) files. [Morpheus](https://software.broadinstitute.org/morpheus/) provides a convenient way to convert your data tables into GCT format.

For more information about the GSEA method and MSigDB please visit http://software.broadinstitute.org/gsea/.

### panoply_ssgsea_report

This module creates an [R Markdown](https://rmarkdown.rstudio.com/) report to provide a high-level summary of the [panoply_ssgsea](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea) module results.

The report provides:

* A heatmap depicting the normalized enrichment scores (NES) of gene sets significant in at least one data column. Significance is defined by parameter `fdr` in section `panoply_ssgsea_report` of the `cfg_yaml` file.

* List of parameters used in the [panoply_ssgsea](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea) module.


## Input

### Required inputs:

* ```input_ds```: (`.gct` file) input GCT file
* ```gene_set_database```: (`.gmt` file) gene set database
* ```yaml_file```: (`.yaml` file) master-parameters.yaml
* ```output_prefix```: (String) File prefix for output files.

* ```preprocess_gct```: (Boolean) If FALSE [panoply_preprocess_gct](./Data-Preparation-Modules%3A-panoply_preprocess_gct) will be skipped and the GCT file will be used as is (default: FALSE).

### Optional inputs:

#### panoply_preprocess_gct
* ```acc_type```: (String) Type of accession number in 'rid' object in GCT file ("uniprot", "refseq" (default), "symbol").
* ```id_type```: (String) Notation of site-ids: 'sm' - Spectrum Mill (default); 'wg' - Web Gestalt; 'ph' - Philosopher. Only relevant for PTM-SEA.
* ```id_type_out``` (String) Type of site id for output: 'uniprot'(default), 'refseq', 'seqwin'. Only relevant for PTM-SEA.
* ```level``` (String) Mode of report:
  + 'ssc' - single-site-centric
  + 'gc' - gene-centric (default)
  + 'gcr' - gene-centric-redundant
* ```loc``` (Boolean) If TRUE only fully localized sites will be considered (default: TRUE). Localization infromation is expected to be encoded in the site identifier. Respective parsing rules are determined by '--id_type'.
* ```gene_col```: (String) Name of column listing gene names; used for gene centric reports (default: "geneSymbol").
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

#### panoply_ssgsea
* ```correl_type```: (String) Correlation type: "z.score" (default), "rank", "symm.rank".
* ```global_fdr```: (Boolean) If TRUE global FDR across all data columns is calculated (default: FALSE).
* ```min_overlap```: (Integer) Minimal overlap between signature and data set (default: 10).
* ```tolerate_min_overlap_err``` (String) true/false toggle for tolerating "not-enough-overlap" errors. Recommended value is FALSE.
* ```nperm```: (Integer) Number of permutations (default: 1000).
* ```output_score_type```: (String) Score type: "ES" - enrichment score,  "NES" - normalized ES (default).
* ```sample_norm_type```: (String) Sample normalization: "rank"(default), "log", "log.rank"
* ```statistic```: (String) Test statistic: "area.under.RES" (default), "Kolmogorov-Smirnov"
* ```weight```: (Float) When weight=0, all genes have the same weight; if weight>0 actual values matter and can change the resulting score (default: 0.75).
* ```output_prefix```: (String, default="results-ssgsea") File prefix for output files.



## Output

`panoply_ssgsea_workflow` produces the follow outputs:

* ```results```: (`.tar`) An output tar file that contains the ssgsea analysis results 
* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
