Documentation at https://github.com/broadinstitute/PANOPLY/blob/version-1_0/release/version-1_0/panoply_cons_clust/panoply_cons_clust.md

# ```panoply_cons_clust```

## Description

This module utilizes resampling-based consensus clustering (Monti et al., 2003) to derive robust proteome clusters. First, the data matrix is filtered to remove all proteins with standard deviation (sd) of <`$clustering_sd_threshold` across sample columns. The resulting data matrix is then transformed into 1,000 bootstrap sample data sets which are clustered into *K* clusters using k-means. From the distribution of bootstrap iterations, a consensus matrix is constructed whose the entries _(i, j)_  record the number of times items _i_ and _j_ were assigned to the same cluster divided by the total number of times both items are selected. Thus, a perfect consensus matrix corresponds to a matrix with all entries equal either 0 or 1. A range of possible cluster numbers *K* between 2 and 10 are evaluated and the best *K* is determined by comparing the *empirical cumulative distribution (CDF)* of the resulting consensus matrices. To compare the clusterings the increase of CDF area *K <sub>delta</sub>* is evaluated and the *K* with the largest *K <sub>delta</sub>* is defined as best *K*.

To identify molecular features separating the *K* clusters, the clustering results are subjected to differential marker selection analysis comparing samples in one cluster to all other samples (*one-vs-other*) using the methods implemented in the [panoply_association](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_association) module. Similariy, Gene Set Enrichment Analysis (GSEA) (Subramanian et al., 2005) on the clustering results is performed (*one-vs-other*).

## Input

### Required inputs:

* ```inputData```: (`.tar` file or `.gct` file) tarball from ```panoply_association``` or other PANOPLY module;\
(`.gct` file) normalized/filtered input if ```standalone``` is ```TRUE```
* ```type```: (String) data type
* ```yaml```: (`.yaml` file) parameters in `yaml` format
* ```standalone```: (String) set to ```TRUE``` to run as a self-contained module;\
if ```TRUE``` the ```analysisDir``` input is required

### Optional inputs:

* ```groupsFile```: (`.csv` file) subset of sample annotations
* ```analysisDir```: (String) name of analysis directory
* ```clustering_sd_threshold```: (Int, default = 2) threshold for filtering data before consensus clustering 
* ```clustering_na_threshold```: (Float, default = 0.5) maximum fraction of missing values for clustering; rest are imputed
* ```outFile```: (String, default = 'panoply_cluster-output.tar') name for output `.tar` file

## Output

* ```outputs```: (File) tarball including:


**Clustering results**

* text files (`.csv`, `.cls`, `.gct`)
  * `${type}-bestclus.csv`: cluster labels for best K
  * `${type}-Cluster.cls`: cluster labels for best K in [cls](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Categorical_.28e.g_tumor_vs_normal.29_class_file_format_.28.2A.cls.29) format
  * `${type}-Cluster-class.X[2-10].cls`: binary cluster labels (other vs. cluster X) for each cluster using best K in [cls](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Categorical_.28e.g_tumor_vs_normal.29_class_file_format_.28.2A.cls.29) format
  * `${type}-cluster-enrichment.csv`: result of an overrepresenatation analysis (Fisher's exact test) of categorial sample annotaions specified in ```groupsFile``` applied to each cluster (best K).
 
* `.pdf` files
  * `${type}_proteome_cluster_metrics_K2-K10.pdf`: cluster metrics (delta AUC, silhouette score, cophenetic correlation) as a function of cluster numbers 
  * `${type}.bestclus.silhouette.pdf`: silhouette plot for best K
  * `${type}_cluster_pca_K[2-10].pdf`: Principal component analyis for best K (PC1 vs PC2)
  * `${type}_consensus_matrix_k[2-10].pdf`: consensus matrix for best K
  
* `.png` files
	* `${type}_consensus_matrix_k[2-10].png`: consensus matrices for cluster numbers 2 to 10

 
**Marker selection results**

* List of significant markers derived using SAM.
  * `${type}-Cluster-analysis-markers-all.csv`: All features from the input data matrix.
  * `${type}-Cluster-analysis-markers-fdr0.01.csv`: All features separating the clusters (SAM FDR<0.01).
  * `${type}-Cluster-analysis-markerimp-fdr0.01.csv`: All features separating the clusters (SAM FDR<0.01) ranked by importance.


* Cluster-specific list of significant markers derived using SAM.
  * `${type}-Cluster-class.X[1-K]-analysis-markers-all.csv`: All features from the input data matrix.
  * `${type}-Cluster-class.X[1-K]-analysis-markers-fdr0.01.csv`: All features separating the clusters (SAM FDR<0.01).
  * `${type}-Cluster-class.X[1-K]-analysis-markerimp-fdr0.01.csv`: All features separating the clusters (SAM FDR<0.01) ranked by importance.


* Summary of model and prediction results
  * `${type}-Cluster-analysis-model-results.txt`: Summary of the different classifiers (Partial Least Squares, PLS; Nearest Shrunken Centroids, NSC; Random Forest, RF; glmnet)
  * `${type}-Cluster-analysis-train-results-fdr0.01.csv`: Table of prediction results for training data from all classifiers.

**GSEA results** 

* `${type}-Cluster-class.X[1-K]-analysis-gsea-analysis` (folder): GSEA outputs of a two class comparison for each cluster (other clusters vs current cluster)

## References

* Monti, S. et al. **Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data.** Mach. Learn. 52, 91-118 (2003).

* Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005).
   **Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.**
  Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. http://doi.org/10.1073/pnas.0506580102

	