# ```panoply_cna_correlation```

## Description

This module runs CNA analysis to measure cis- and trans-correlations of CNA with RNA and protein/PTM expression levels using the [WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html) correlation calculation library. Protein expression data is mapped to gene symbols, and harmonized data tables for CNA, RNA and protein/PTM are created, with all data types having the same set of (common) genes and samples (in the `panoply_cna_setup` module). Cis-correlation calculates CNA-RNA and CNA-protein correlation for each gene; trans-correlations determine CNA-RNA and CNA-protein correlations for *all pairs of genes*. P-values are also calculated for each correlation. Significant CNA-RNA and CNA-protein (positive or negative) correlations are identified after correcting for multiple testing, and plotted on CNA-vs-RNA and CNA-vs-protein correlations plots with genes organized along the axes by chromosome location. See (Zhang et al., 2014) and (Mertins et al., 2016).

## Input

Require inputs:

* ```tarball```: (`.tar` file) tarball from ```panoply_cna_setup```, containing harmonized data tables for CNA, RNA and protein/PTM
* ```type```: (String) proteome data type 
* ```yaml```: (`.yaml` file) parameters file in `yaml` format

Optional inputs:

* ```fdr_cna_corr```: (Float) FDR p-value for determining significant cis- and trans-correlations
* ```outFile```: (String, default = "panoply_cna_correlation-output.tar") output `.tar` file name


## Output

Tarball including the following files in the `cna` subdirectory:

* CNA-RNA correlation:
	* `all-mrna-vs-cna-corr.csv` gene x gene table of CNA-RNA correlation values 
	* `all-mrna-vs-cna-pval.csv` table of CNA-RNA multiple testing corrected FDR p-values for corresponding correlation values
	* `all-mrna-vs-cna-sigevents.csv` count of significant correlations ("events") for each CNA gene; the table also includes an indication of whether the CNA-RNA cis-correlation was statistically significant, and the chromosome location for each gene
* CNA-proteome correlation:
	* `all-pome-vs-cna-corr.csv` gene x gene table of CNA-proteome (protein or PTM) correlation values 
	* `all-pome-vs-cna-pval.csv` table of CNA-RNA multiple testing corrected FDR p-values for corresponding correlation values
	* `all-pome-vs-cna-sigevents.csv` count of significant correlations ("events") for each CNA gene; the table also includes an indication of whether the CNA-proteome (protein or PTM)q cis-correlation was statistically significant, and the chromosome location for each gene

* CNA correlation plot (`all-cna-plot.png`) plot showing correlation between (i) CNA and RNA expression and between (ii) CNA and protein/PTM abundance. Significant positive and negative correlations are indicated in red and green, respectively. CNA-driven cis effects usually appear as a red diagonal line; trans effects appear as vertical lines with red/green dots. The accompanying histograms at the bottom of each correlation plot show the number of significant cis and trans events corresponding to the indicated genomic loci (upward plot) as well as the overlap between CNA-RNA and CNA-protein events (downward plot).

## References

* Zhang, B., Wang, J., XiaojingWang, Zhu, J., Liu, Q., Shi, Z., Chambers, M., Zimmerman, L., Shaddox, K., Kim, S.,et al. (2014). Proteogenomic characterization of human colon and rectal cancer. *Nature* https://dx.doi.org/10.1038/nature13438
* Mertins, P., Mani, D., Ruggles, K., Gillette, M., Clauser, K., Wang, P., Wang, X., Qiao, J., Cao, S., Petralia, F., et al. (2016). Proteogenomics connects somatic mutations to signalling in breast cancer. *Nature*  534(7605), 55 - 62. https://dx.doi.org/10.1038/nature18003.

