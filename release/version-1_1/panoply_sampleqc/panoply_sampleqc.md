# ```panoply_sampleqc```

## Description

This module assess sample QC using CNA, RNA and proteome data (mapped to gene symbols) based on the following methods (Mertins et al., 2016):

* Sample-level correlation between proteome, RNA, and CNA. Inspecting heatmaps of pairwise correlation between these data types can identify sample swaps or mislabeling
* Co-clustering of RNA and protein data, using only genes with good correlation (> `corThreshold`). Similarity of sample material used for RNA and protein profiling should result in RNA and protein clustering together for a large fraction of samples
* Comparison of ESTIMATE (Yoshihara et al., 2013) scores for RNA, protein and CNA.

## Input

Required inputs:

* ```tarball```: (`.tar` file) tarball from ```panoply_harmonize```
* ```type```: (String) proteomics data type
* ```yaml```: (`.yaml` file) parameters in `yaml` format

Optional inputs:
* ```corThreshold```: (Float, default = 0.4) correlation threshold for filtering genes in RNA and protein data for co-clustering
* ```outFile```: (String, default = "panoply_sampleqc-output.tar") output `.tar` file name

## Output

Tarball including the following files in the `sample-qc` subdirectory: 

* ESTIMATE score table for RNA (`rna-estimate-scores.gct`), CNA (`cna-estimate-scores.gct`) and proteome (`*-estimate-scores.gct`)
* Plots showing correlation heatmaps, co-clustering fanplot and boxplots for ESTIMATE scores (`sample-qc-plots.pdf`)
	
## References

* Mertins, P., Mani, D., Ruggles, K., Gillette, M., Clauser, K., Wang, P., Wang, X., Qiao, J., Cao, S., Petralia, F., et al. (2016). Proteogenomics connects somatic mutations to signalling in breast cancer. *Nature*  534(7605), 55 - 62. https://dx.doi.org/10.1038/nature18003.
* Yoshihara, K., Shahmoradgoli, M., nez, E., Vegesna, R., Kim, H., Torres-Garcia, W., o, V., Shen, H., Laird, P., Levine, D., Carter, S., Getz, G., Stemke-Hale, K., Mills, G., Verhaak, R. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. *Nature Communications*  4(), 1 - 11. https://dx.doi.org/10.1038/ncomms3612.
