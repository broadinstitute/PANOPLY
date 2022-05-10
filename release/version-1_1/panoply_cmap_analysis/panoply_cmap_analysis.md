# ```panoply_cmap_analysis```

## Description
CMAP analysis identifies candidate genes driving response to copy number alternations using large-scale Connectibity Map ([CMAP](https://clue.io/cmap)).

The CMAP (Lamb et al., 2006; Subramanian et al., 2017) is a collection of about 1.3 million gene expression profiles from cell lines treated with bioactive small molecules (~20,000 drug perturbagens), shRNA gene knockdowns (~4,300) and ectopic expression of genes. The CMAP dataset is available on GEO ([Series GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)). For this analysis, we use the Level 5 (signatures from aggregating replicates) TouchStone dataset with 473,647 total profiles, containing 36,720 gene knock-down profiles, with measurements for 12,328 genes. See [CLUE](https://clue.io/GEO-guide) for more information. This dataset is split into 100 fragments and stored in a [Google bucket](gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data) in order to allow parallel processing of queries.

To identify candidate driver genes, proteome profiles of copy number-altered samples are correlated with gene knockdown mRNA profiles in the above CMAP dataset, and enrichment of up/downregulated genes is evaluated. Normalized log2 copy number values less than -`cna_threshold` define deletion (loss), and values greater than +`cna_threshold` define copy number amplifications (gains). In the copy number-altered samples (separately for CNA amplification and CNA deletion), the trans-genes (identified by significant correlation in `panoply_cna_correlation` with pvalues less than `fdr_pvalue`) are grouped into UP and DOWN categories by comparing the protein ratios of these genes to their ratios in the copy number neutral samples (normalized log2 copy number between -`cna_threshold` and +`cna_threshold`). The lists of UP and DOWN trans-genes are then used as queries to interrogate CMAP signatures and calculate weighted connectivity scores (WTCS) using the single-sample GSEA algorithm in `panoply_ssgsea`. The weighted connectivity scores are then normalized for each perturbation type and cell line to obtain normalized connectivity scores (NCS). See (Subramanian et al., 2017) for details on WTCS and NCS. For each query we then identify outlier NCS scores, where a score is considered an outlier if it is beyond 1.5 times the interquartile range of score (IQR) distribution for the query. The query gene is designated a candidate driver if (i) the score outliers are statistically cis-enriched (Fisher test with BH-FDR multiple testing correction p-value less than `cmap_fdr`) and (ii) the gene has statistically significant and positive cis-correlation. **NOTE** When `legacy_score` is set to `TRUE`, an older scoring method described in (Mertins et al., 2016) is used.

For a gene to be considered for inclusion in a CMAP query it needs to i) have a copy number change (amplification or deletion) in at least `cna_effects_threshold` samples; ii) have at least `min_sigevents` significant trans genes; and iii) be on the list of shRNA knockdowns in the CMAP. CMAP queries on then run on the `top.N` genes with the most trans-events. 

In order to ensure that the identified candidate driver genes are not a random occurrence, a permutation test (with `n_permutations`) is used to determine how many candidate driver genes would be identified with random input (Gillette et al., 2020). To determine FDR, each permutation run is treated as a Poisson sample with rate l, counting the number of identified candidate driver genes. A Score confidence interval is calculated and the midpoint of the confidence interval used to estimate the expected number of false positives. 

To identify how many trans-correlated genes for all candidate regulatory genes can be directly explained by gene expression changes measured in the CMAP shRNA perturbation experiments, knockdown gene expression consensus signature z-scores (knockdown/control) are used to identify regulated genes, followed by counting the number of trans-genes in this list of regulated genes.

To obtain biological insight into the list of candidate driver genes, we perform (i) enrichment analysis on samples with extreme CNA values (amplification or deletion) to identify statistically enriched sample annotation subgroups from annotation provided in the `cmap_enrichment_groups` file; and (ii) GSEA on cis/trans-correlation values to find enriched pathways in `annotation_pathway_db`.

## Input

Required inputs:

* ```CNAcorr_tarball```: (`.tar` file) tar file from output of `panoply_cna_correlation` that contains harmonized RNA, CNA and proteomics data, along with CNA cis/trans correlation tables with corresponding p-values.
* ```subset_list_file```: (File) list of CMAP dataset fragments, default set to gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/cmap-data-subsets-index.txt
* ```cmap_level5_data```: (File) complete CMAP dataset, default location gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/annotated_GSE92742_Broad_LINCS_Level5_COMPZ_geneKDsubset_n36720x12328.gctx
* ```subset_bucket```: (String) location of subdirectory containing CMAP dataset fragments, default set to gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/cmap-data-subsets
* ```annotation_pathway_db```: (File) pathway database for GSEA enrichment analysis to promote biological insight into list of candidate driver genes
* ```n_permutations```: (Int) number of permutations to run for determining overall FDR 
* ```yaml```: (File) parameters file in `yaml` format

Optional inputs:

* ```group```: (String, default = "all") filename prefix used for files created during CMAP analysis
* ```data_type```: (String, default = "pome") omics data type; suppoted types are "pome" (all proteomics/PTM types) and "mrna" (for RNAseq data)
* ```cmap_enrichment_groups```: (`.csv` File) subset of sample annotations, providing classes for enrichment analysis of candidate genes
* ```cna_threshold```: (Float, default = 0.3) copy number up/down threshold; copy number is considered UP regulated if > `cna_threshold` and DOWN regulated if < `-cna_threshold` 
* ```cna_effects_threshold```: (Int, default = 15) minimum number of samples with up/down copy number that must be present to include gene for CMAP analysis
* ```min_sigevents```: (Int, default = 20) gene must have at least this many significant trans events to be considered for inclusion
* ```max_sigevents```: (Int, default = 1800) if a gene has more then `max.sigevents` trans events, the top `max.sigevents` will be chosen
* ```top_N```: (Int, default = 500) maximum number of genes to run CMAP analysis on
* ```fdr_pvalue```: (Float, devault = 0.05) FDR for CNA correlation significance level
* ```log_transform```: (String, default = "FALSE") if TRUE, log2 transform input data
* ```must_include_genes```: (String) an R vector in `c ('gene1', 'gene2', ...)` notation to list genes that must be included in the CMAP analysis
* ```cis_fdr```: (String, default = fdr_pvalue) FDR for cis-correlation significance level 
* ```legacy_score```: (String, default = "FALSE") if TRUE, legacy connectivity score will be calculated (using mean rank points), with permutation FDR (see Description)
* ```rankpt_n```: (String, default = 4) number of CMAP profiles to consider when calculating mean rank point (used only when `legacy_score` is `TRUE`)
* ```mean_rankpt_threshold```: (Int, default = 85) minimum value of mean rank point for gene signature to be considered enriched (used only when `legacy_score` is `TRUE`)
* ```cmap_fdr```: (Float, default = 0.25) FDR threshold for Fisher test on outlier scores in order for gene to be considered a candidate driver
* ```alpha```: (String, default = 0.05) p-value threshold for CMAP profile z-scores (see Description) and enrichments

**NOTE:** The inputs described above are for the primary workflow. Additional optional inputs for tasks constituting the workflow are already set to appropriate defaults and do not need to be modified.


## Output

Tarball including the following files in the `cmap` subdirectory:

* A list of significant candidate driver genes identified using the Fisher test on outlier scores for 
  - genes with CNA UP **or** DOWN (unidirectional, `*-sig-genes-unidirectional.txt`)
  - genes with CNA UP **and** DOWN (bidirectional, `*-sig-genes-bidirectional.txt`)
* FDR score for uni- and bidirectional significant candidate driver genes (`*-sig-genes-with-fdr.txt`)
* A gene x sample table indicating outlier status (`*-outliers.csv`)
* List and plot of extreme samples (with CNA UP/DOWN) for each significant candidate driver gene (`*-sig-genes-extreme-samples.*`)
* Enrichment these extreme samples in various sample annotation groups (`*-sig-genes-enrichment.csv` and `*-sig-genes-enrichment-pval*.csv` with the latter showing only group that are enriched with `p.value < alpha`
* Overlap of trans genes with extreme genes in CMAP profiles (based on z-score); `*-sig-genes-overlap.gmt` lists overlap for each significant candidate driver gene and these results are graphically shown in `*-sig-genes-overlap.pdf`
* Input (`*-permuted-genes-???.gmt`) and output results (`*-permutation-*`) for all permutations 

## References
* Lamb, J., Crawford, E.D., Peck, D., Modell, J.W., Blat, I.C., Wrobel, M.J., Lerner, J., Brunet, J.-P., Subramanian, A., Ross, K.N., et al. (2006). The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. *Science* 313, 1929–1935.
* Subramanian, A., Narayan, R., Corsello, S.M., Peck, D.D., Natoli, T.E., Lu, X., Gould, J., Davis, J.F., Tubelli, A.A., Asiedu, J.K., et al. (2017). A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell* 171, 1437–1452.e17.
* Gillette, M., Satpathy, S., Cao, S., Dhanasekaran, S., Vasaikar, S., Krug, K., Petralia, F., Li, Y., Liang, W., Reva, B., et. al. (2020). Proteogenomic Characterization Reveals Therapeutic Vulnerabilities in Lung Adenocarcinoma. *Cel*l  182(1), 200 - 225.e35.
* Mertins, P., Mani, D., Ruggles, K., Gillette, M., Clauser, K., Wang, P., Wang, X., Qiao, J., Cao, S., Petralia, F., et al. (2016). Proteogenomics connects somatic mutations to signalling in breast cancer. *Nature*  534(7605), 55 - 62. 
