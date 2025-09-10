# ```panoply_metaboanalyst```

## Description

This module is a PANOPLY adaptation of [MetaboAnalyst](https://www.metaboanalyst.ca/MetaboAnalyst/). Specifically, this module adapts the Joint Pathway-Analysis tool, first described under Integrated Pathway Analysis in [Xia, J. et al. 2015](https://doi.org/10.1093/nar/gkv380). The tool takes both a metabolomic dataset and a gene-level dataset (proteomic or transcriptomic), and performs combined pathway enrichment across a curated list of KEGG pathways (pulled 12/20/2023 via KEGG API; (Pang, Z. et al 2024)). In order to combine pathway enrichment results from two datasets, overreresentation analysis is performed separately on metabolites and genes, and p-values are combined using one of four methods (`pval_comb`); additional documentation on these p-value combination methods can be found on [the MetaboAnalyst webpage](https://www.metaboanalyst.ca/MetaboAnalyst/). In addition to performing overrepresentation analysis, topological analysis is performed to calculate the impact of the significantly altered features on each pathway, and pathway-visualizations are created to show significant features.

The PANOPLY version of MetaboAnalyst extends this tool by automating analysis across multiple annotations of interest. The user should provide a metabolomic and gene-level (proteomic or transcriptomic) dataset, as well as a `groups_file` with all sample-annotations to be analyzed. For each annotation value, a 2-class (one-vs-rest) mod-T-test will be performed on metabolites and gene-level features, to get a list of significantly altered features; if there are an appropriate number of significant features, overrepresentation analysis and topological analysis are performed. This analysis is repeated across all annotations in the provided `groups_file`.

### Metabolite-Only Analysis:

While the main function of this module is to perform combined analysis of proteomic and metabolomic data, support has also been added for MetaboAnalyst's metabolite-only enrichment analysis. If only a metabolomic dataset is provided, pathway enrichment analysis will instead be performed on a curated set of pathways from the Small Molecular Pathway Database (Xia, J. & Wishart D. S. 2010).

The module also supports quantitative enrichment analysis (`anal_type="QEA"`) as an alternative to overrepresentation analysis. In QEA, a global test is performed across all features in a pathway, which provides increased statistical power; the analysis is described in full in [Xia, J. & Wishart D. S. 2010](https://doi.org/10.1093/nar/gkq329). We recommend using QEA _only_ if running metabolomics alone, as it has not been extensively tested for joint-pathway analysis. Topological analysis will not be performed when running quantitative-enrichment analysis, as individual feature-significance is not calculated.

If you anticipate running [panoply_metaboanalyst_report](./Report-Modules%3A-panoply_metaboanalyst_report.md) to generate an interactive report, it is recommended to run Joint-Pathway Analysis by providing both metabolite and an additional -omic dataset, and choosing `anal_type="ORA"`.

## Input

### Required inputs:

* ```meta_gct```: (`.gct` file) normalized metabolomic data table in the form of a GCT
* ```omic_gct```: (`.gct` file) normalized -omic data table in the form of a GCT (**optional** but recommended)
* ```ome_type```: (String) type of -omic data provided, e.g. "RNA" or "proteome" (**optional** but recommended)
* ```groups_file```: (`.csv` file) annotation file, subsetted to annotations of interest for this analysis (**optional** but recommended)

* ```output_prefix```: (String, default="results_metaboanalyst") prefix used to name the output tar file
* ```yaml_file```: (`.yaml` file) master-parameters.yaml

### Optional inputs:

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

* ```results```: `.tar` file containing all results from MetaboAnalyst analysis 
	* For each annotation-of-interest:
		* `.xsl` file containing all pathway-enrichment results for each annotation subvalue
		* Heatmap (`.png` and `.pdf`) summarizing the top significant pathways across all annotation-subvalues
		* `results_<annotation>/` folder containing a top-hits summary figure and network visualizations for the top significant-pathways, for each annotation subvalue
	* **Log File:** (`${}_with${type}_log_file.csv`) summary table of all annotations analyzed; lists whether an annotation value was valid (i.e. is not NA or blank) and whether it had valid results (i.e. significant features and significant pathway-enrichments).


## References

1. Pang, Z. et al. MetaboAnalyst 6.0: towards a unified platform for metabolomics data processing, analysis and interpretation. Nucleic Acids Research 52, W398–W406 (2024).

2. Xia, J., Sinelnikov, I. V., Han, B. & Wishart, D. S. MetaboAnalyst 3.0—making metabolomics more meaningful. Nucleic Acids Res 43, W251–W257 (2015).

3. Xia, J. & Wishart, D. S. MSEA: a web-based tool to identify biologically meaningful patterns in quantitative metabolomic data. Nucleic Acids Research 38, W71–W77 (2010).

