# ```panoply_metaboanalyst_workflow```

## Description

This workflow executes a PANOPLY adaptation of [MetaboAnalyst](https://www.metaboanalyst.ca/MetaboAnalyst/), performing [pathway enrichment](./Data-Analysis-Modules%3A-panoply_metaboanalyst) of metabolomic (and proteomic/transcriptomic) features and generates an [interactive report](./Report-Modules%3A-panoply_metaboanalyst_report.md). The workflow is based off the 6.0 release of MetaboAnalyst; all databases are pulled from that version (Pang, Z. et al. 2024). It executes the following modules:

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_metaboanalyst</code>](./Data-Analysis-Modules%3A-panoply_metaboanalyst)         |  performs joint-pathway enrichment analysis on metabolites and an additional -ome  |
| [<code>panoply_metaboanalyst_report</code>](./Report-Modules%3A-panoply_metaboanalyst_report.md)   |  creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the pathway-analysis results |


## Input

### Required inputs:

* ```meta_gct```: (`.gct` file) normalized metabolomic data table in the form of a GCT
* ```omic_gct```: (`.gct` file) normalized -omic data table in the form of a GCT (**optional** but recommended)
* ```ome_type```: (String) type of -omic data provided, e.g. "RNA" or "proteome" (**optional** but recommended)
* ```groups_file```: (`.csv` file) annotation file, subsetted to annotations of interest for this analysis (**optional** but recommended)

* ```output_prefix```: (String) prefix used to name the results files
* ```yaml_file```: (`.yaml` file) master-parameters.yaml

### Optional inputs:

#### panoply_metaboanalyst

* ```anal_type```: (String, default="ORA") type of pathway-enrichment analysis to be performed ("ORA" for Overrepresentation Analysis or "QEA" for Quantitative Enrichment Analysis); **ORA is recommended** for joint-pathway analysis
* ```pval_comb```: (String, default="pvalo") method for combining p-values in multiomic enrichment analysis; options include "query" (combine queries), "pvalu" (unweighted), "pvalo" (overall), or "pvalp" (pathway-level)

* ```max_annot_levels```: (Int, default=10) max number of levels in an annotation category, to be considered qualitative and included in the analysis
* ```pval_signif```: (Float, default=0.05) p-value threshold for significant enrichement
* ```top_n_networks```: (Int, default=10) top N networks to plot per annot subvalue

* ```meta_id_col```: (String, default="metab_ids") rdesc column in Metabolite GCT with meta_id_type IDs; if set to `NULL`, GCT rid values will be used.
* ```meta_id_type```: (String, default="hmdb_id") ID type used for metabolites
* ```gene_column```: (String, default="geneSymbol")  rdesc column in -omic GCT with gene IDs
* ```gene_id_type```: (String, default="SYMBOL") ID type used for gene IDs


## Output

`panoply_metaboanalyst_workflow` produces the follow outputs:

* ```metaboanalyst_tar```: (`.tar`) An output tar file that contains the metaboanalyst analysis results 
* ```metaboanalyst_report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.


## References

1. Pang, Z. et al. MetaboAnalyst 6.0: towards a unified platform for metabolomics data processing, analysis and interpretation. Nucleic Acids Research 52, W398–W406 (2024).
