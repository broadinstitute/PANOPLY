Documentation at https://github.com/broadinstitute/PANOPLY/blob/release-1_0/release/version-1_0/panoply_rna_protein_correlation/panoply_rna_protein_correlation.md

# ```panoply_rna_protein_correlation```

## Description

Measures and plots correlation between mRNA expression and protein abundance for each gene-protein pair using [pearson](https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/cor.test) correlation. Protein IDs are mapped to gene symbols; PTM sites are rolled-up to the protein level before mapping to gene symbols. The module outputs correlation tables and histograms summarizing RNA-protein correlation. See (Zhang et al., 2014) and (Mertins et al., 2016).

## Input

Required inputs:

* ```inputData```: (`.tar` file) tarball from ```panoply_normalize_ms_data``` or normalized protein/PTM site data in `gct` format (when `standalone` is `TRUE`)
* ```rnaExpr```: (`.gct` file) RNA expression data
* ```type```: proteomics data type
* ```standalone```: (String) set to ```TRUE``` to run as a self-contained module; if ```TRUE``` the ```analysisDir``` input is required
* ```yaml```: (`.yaml` file) parameters in `yaml` format

Optional inputs:

* ```rnaSDthreshold```: (Int, default = 1) for standard deviation variation filter; set to NA to disable
* ```profilePlotTopN```: (Int, default = 25) plots showing RNA and protein levels across samples for the `topN` gene-protein pairs (sorted by correlation) are generated
* ```analysisDir```: (String) name of analysis directory
* ```outFile```: (String, default = "panoply_rna_protein_correlation-output.tar") output `.tar` file name

## Output

Tarball including the following files in the `rna` subdirectory: 

* Harmonized (`rna-seq.gct`) and filtered (`rna-seq-sdfilter.gct`) RNA data with sample order matching that in the proteome data table
* Tables listing RNA-protein correlation for every gene-protein pair (`*-mrna-cor.tsv`) and the best gene-protein pairs with highest (correlation > 0.7), statistically significant correlation (`*-mrna-cor-best.tsv`)
* Plots showing histograms of RNA-protein correlation for the following gene-protein pairs:
	* all gene-protein pairs (`*-mrna-cor.pdf`)
	* best pairs (`*-mrna-cor-best.pdf`)
	* statistically significant pairs (`*-mrna-cor-sig.pdf`)
	* all pairs, with statistically significant pairs highlighted (`*-mrna-cor-combined.pdf`)
	* plots showing RNA and protein levels across samples for the `topN` gene-protein pairs

## References

* Zhang, B., Wang, J., XiaojingWang, Zhu, J., Liu, Q., Shi, Z., Chambers, M., Zimmerman, L., Shaddox, K., Kim, S.,et al. (2014). Proteogenomic characterization of human colon and rectal cancer. *Nature* https://dx.doi.org/10.1038/nature13438
* Mertins, P., Mani, D., Ruggles, K., Gillette, M., Clauser, K., Wang, P., Wang, X., Qiao, J., Cao, S., Petralia, F., et al. (2016). Proteogenomics connects somatic mutations to signalling in breast cancer. *Nature*  534(7605), 55 - 62. https://dx.doi.org/10.1038/nature18003.
