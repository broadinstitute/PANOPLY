# PANOPLY: A cloud-based platform for automated and reproducible proteogenomic data analysis
#### Version 1.0

PANOPLY is a platform for applying state-of-the-art statistical and machine learning algorithms to transform multi-omic data from cancer samples into biologically meaningful and interpretable results. PANOPLY leverages [Terra](http://app.terra.bio)—a cloud-native platform for extreme-scale data analysis, sharing, and collaboration—to host proteogenomic workflows, and is designed to be flexible, automated, reproducible, scalable, and secure. A wide array of algorithms applicable to all cancer types have been implemented, and available in PANOPLY for analysis of cancer proteogenomic data.


![*Figure 1.* PANOPLY architecture overview. Inputs to PANOPLY consists of (i) externally characterized genomics and proteomics data (in gct format); (ii) sample phenotype and annotations (in csv format); and (iii) parameters settings (in yaml format). Input data are subject to quality checks followed by optional normalization and filtering for proteome data. Analysis ready data tables are then used as inputs to the data analysis modules. Results from the data analysis modules are summarized in interactive reports.](panoply-overview.png)


PANOPLY v1.0 consists of the following components:

* A Terra production workspace on [PANOPLY_Production_Pipelines_v1_0](https://app.terra.bio/#workspaces/broad-firecloud-cptac/PANOPLY_Production_Pipelines_v1_0) with a preconfigured unified workflow to automatically run all analysis tasks on proteomics (global proteome, phosphoproteome, acetylome, ubiquitylome), transcriptome and copy number data. Figure 1 provides an overview of PANOPLY architecture, and the various tasks that constitute the complete workflow.
* A second workspace [PANOPLY_Production_v1_0](https://app.terra.bio/#workspaces/broad-firecloud-cptac/PANOPLY_Production_v1_0) that includes separate methods for each analysis component of the unified PANOPLY pipeline. These methods can be run independently, or combined into custom workflows using WDL scripts.
* An interactive notebook included with the production workspaces that provides step-by-step instructions for uploading data, identifying data types, specifying parameters, and setting up the PANOPLY workspace.
* A GitHub [repository](https://github.com/broadinstitute/PANOPLY) that contains code for all PANOPLY tasks and workflows, including R code for implementing analysis algorithms, task module creation, and release management.
* A GitHub [wiki](https://github.com/broadinstitute/PANOPLY/wiki) with documentation and a tutorial.

PANOPLY provides a comprehensive collection of proteogenomic data analysis methods including sample QC (sample quality evaluation using profile plots and tumor purity scores<sup>1</sup>, identify sample swaps, etc.), association analysis, RNA and copy number correlation (to proteome), connectivity map analysis<sup>1 ,2</sup>, outlier analysis using BlackSheep<sup>3</sup>, PTM-SEA<sup>4</sup>, GSEA<sup>5</sup> and single-sample GSEA<sup>6</sup>, consensus clustering, and multi-omic clustering using non-negative matrix factorization (NMF). Most analysis modules include a report generation task that outputs a HTML interactive report summarizing results from the respective analysis tasks. Mass spectrometry-based proteomics data amenable to analysis by PANOPLY includes isobaric label-based LC-MS/MS approaches like iTRAQ, TMT and TMTPro as well as label-free proteomics data at proteome and multiple PTM-omes including phospho-, acetyl-and ubiquitylomes.

Users can also quickly [add their own tasks](https://support.terra.bio/hc/en-us/articles/360031366091-Create-edit-and-share-a-new-workflow) and integrate them into the PANOPLY workflow.


## Quick Start
For a quick introduction and tour of PANOPLY, follow the [tutorial](https://github.com/broadinstitute/PANOPLY/wiki/PANOPLY-Tutorial). Detailed documentation can be found [here](https://github.com/broadinstitute/PANOPLY/wiki).

### Contact

Email proteogenomics@broadinstitute.org with questions, comments or feedback.


### References

1. Mertins, P. et al. Proteogenomics connects somatic mutations to signalling in breast cancer. Nature 534, 55–62 (2016).
2. Subramanian, A. et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell 171, 1437–1452.e17 (2017).
3. Blumenberg, L. et al. BlackSheep: A Bioconductor and Bioconda package for differential extreme value analysis. doi:10.1101/825067.
4.	Krug, K. et al. A Curated Resource for Phosphosite-specific Signature Analysis. Mol. Cell. Proteomics 18, 576–593 (2019).
5.	Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. U. S. A. 102, 15545–15550 (2005).
6.	Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y. & Moody, S. E. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature (2009).
