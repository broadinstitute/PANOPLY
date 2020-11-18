Documentation at https://github.com/broadinstitute/PANOPLY/blob/version-1_0/release/version-1_0/panoply_immune_analysis/panoply_immune_analysis.md

# ```panoply_immune_analysis```

## Description

This module runs several algorithms assigning immune scores for understanding the tumor microenvironment, based on RNAseq data. **E**stimation of **ST**romal and **I**mmune cells in **MA**lignant **T**umor tissues using **E**xpression data (**ESTIMATE**, Yoshihara et al., 2013) uses gene expression signatures to calculate the fraction of stromal and immune cells and infer tumor purity. **xCell** (Aran et al., 2017) uses a gene signature-based method to infer immune and stromal cell types. **ImmuneSubtypeClassifier** (Thorrson et al., 2018) uses immune gene expression signatures to classify tumor samples into one of 6 immune subtypes. Enrichment analysis (Fisher\'s exact test) is performed on immune subtypes to determine enrichment of sample annotation groups specified in the `groupsFile`.


## Input

Required inputs:

* ```inputData```: (File) tarball output (`.tar`) for one of the other modules, or RNAseq data in `gct` format (when `standalone` is `TRUE`)
* ```type```: (String) proteome data type
* ```standalone```: (String) set to ```TRUE``` to run as a self-contained module; if ```TRUE``` the ```analysisDir``` input is required
* ```yaml```: (`.yaml` file) parameters in `yaml` format

Optional inputs:

* ```analysisDir```: (String) name of analysis directory
* ```groupsFile```: (`.csv` file) subset of sample annotations, providing classes for enrichment analysis
* ```fdr```: (Float, default = 0.05) FDR cutoff value for significance in enrichment analysis
* ```heatmapWidth```: (Int, default = 10) width of heatmap generated 
* ```heatmapHeight```: (Int, default = 15) height of heatmap generated
* ```outFile```: (String, default = "panoply_immune_analysis-output.tar") output `.tar` file name



## Output

Tarball including the following within the `immune-analysis` subdirectory:

* `.csv` files:
	* ESTIMATE results (`estimate-scores.csv`)
	* xCell results (`xcell-scores.csv`)
	* ImmuneSubtypeClassifier results (`immune-subtype.csv`)
	* Immune subtype enrichment analysis results (`immune-subtype-enrichment.csv`; filtered results for ```fdr``` in `immune-subtype-enrichment-pval*.csv`)
	
* `.pdf` files:
	* Heatmap of xCell scores (`xcell-scores-heatmap.pdf`)
	* Plots of xCell and ESTIMATE scores (`xCell-vs-ESTIMATE-plots.pdf`)

## References

* Yoshihara, K., Shahmoradgoli, M., nez, E., Vegesna, R., Kim, H., Torres-Garcia, W., o, V., Shen, H., Laird, P., Levine, D., Carter, S., Getz, G., Stemke-Hale, K., Mills, G., Verhaak, R. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. *Nature Communications*  4(), 1 - 11. https://dx.doi.org/10.1038/ncomms3612.
* Aran, D., Hu, Z., Butte, A. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. *Genome Biology* 18(1), 220. https://dx.doi.org/10.1186/s13059-017-1349-1.
* Thorsson, V., Gibbs, D., Brown, S., Wolf, D., Bortone, D., Yang, T., Porta-Pardo, E., Gao, G., Plaisier, C., Eddy, J., et al. (2018). The Immune Landscape of Cancer. *Immunity*  48(4), 812-830.e14. https://dx.doi.org/10.1016/j.immuni.2018.03.023.
