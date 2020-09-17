# PANOPLY: A cloud-based platform for automated and reproducible proteogenomic data analysis
#### Version 1.0

PANOPLY is a platform for applying state-of-the-art statistical and machine learning algorithms to transform multi-omic data from cancer samples into biologically meaningful and interpretable results. PANOPLY leverages [Terra](http://app.terra.bio)—a cloud-native platform for extreme-scale data analysis, sharing, and collaboration—to host proteogenomic workflows, and is designed to be flexible, automated, reproducible, scalable, and secure. A wide array of algorithms applicable to all cancer types have been implemented, and we highlight the application of PANOPLY to the analysis of cancer proteogenomic data.

Using the underlying architecture of Terra, PANOPLY analysis algorithms and methods are implemented as *tasks*, interchangeably referred to as *modules*. A task consists of a command script that is executed in a specified [Docker](https://www.docker.com) container with and that has predefined inputs and outputs. A series of tasks constitutes a *workflow*, that represents a pipeline with outputs of one or more tasks feeding into inputs of downstream tasks. Tasks and workflows are defined using the Workflow Description Language ([WDL](https://support.terra.bio/hc/en-us/sections/360007274612-WDL-Documentation)). 

The central element of PANOPLY is the [Terra workspace](https://support.terra.bio/hc/en-us/categories/360002177552). A workspace is a shareable collection of everything needed for a data analysis project, including data located on cloud storage linked with the workspace, workflows encapsulating algorithms and pipelines, analysis parameters, results, and a data model that organizes sample meta-data into collections of participants, samples, and sample sets. With versioned workflows, workspaces enable reproducible digital research by completely encapsulating input data, analysis methods, settings, and results in a shareable cloud-based space.

![*Figure 1.* PANOPLY architecture overview. Inputs to PANOPLY consists of (i) externally characterized genomics and proteomics data (in gct format); (ii) sample phenotype and annotations (in csv format); and (iii) parameters settings (in yaml format). Input data are subject to quality checks followed by optional normalization and filtering for proteome data. Analysis ready data tables are then used as inputs to the data analysis modules. Results from the data analysis modules are summarized in interactive reports.](panoply-overview.png)


PANOPLY v1.0 consists of the following components:

* A production workspace [PANOPLY_Production_Pipelines_v1_0](https://app.terra.bio/#workspaces/broad-firecloud-cptac/PANOPLY_Production_Pipelines_v1_0) with a preconfigured unified workflow to automatically run all analysis tasks on proteomics (global proteome, phosphoproteome, acetylome, ubiquitylome), transcriptome and copy number data. Figure 1 provides an overview of PANOPLY architecture, and the various tasks that constitute the complete workflow.
* A second workspace [PANOPLY_Production_v1_0](https://app.terra.bio/#workspaces/broad-firecloud-cptac/PANOPLY_Production_v1_0) that includes separate methods for each analysis component of the unified PANOPLY pipeline. These methods can be run independently, or combined into custom workflows using WDL scripts.
* An interactive notebook included with the production workspaces that provides step-by-step instructions for uploading data, identifying data types, specifying parameters, and setting up the PANOPLY workspace.
* A GitHub [repository](https://github.com/broadinstitute/PANOPLY) that contains code for all PANOPLY tasks and workflows, including R code for implementing analysis algorithms, task module creation, and release management.
* A GitHub [wiki](https://github.com/broadinstitute/PANOPLY/wiki) with documentation and a tutorial.

PANOPLY provides a comprehensive collection of proteogenomic data analysis methods including sample QC (sample quality evaluation using profile plots and tumor purity scores3, identify sample swaps, etc.), association analysis, RNA and copy number correlation (to proteome), connectivity map analysis3,26, outlier analysis using BlackSheep27, PTM-SEA28,GSEA29 and single-sample GSEA30, consensus clustering, and multi-omic clustering using non-negative matrix factorization (NMF). Most analysis modules include a report generation task that outputs a HTML interactive report summarizing results from the respective analysis tasks. A complete list of PANOPLY task modules with documentation can be found at the PANOPLY wiki20. Mass spectrometry-based proteomics data amenable to analysis by PANOPLY includes isobaric label-based LC-MS/MS approaches like iTRAQ, TMT and TMTPro as well as label-free proteomics data at proteome and multiple PTM-omes including phospho-, acetyl-and ubiquitylomes.

Users can also quickly [add their own tasks](https://support.terra.bio/hc/en-us/articles/360031366091-Create-edit-and-share-a-new-workflow) and integrate them into the PANOPLY workflow.


## Quick Start
For a quick introduction and tour of PANOPLY, follow the [tutorial](https://github.com/broadinstitute/PANOPLY/wiki/PANOPLY-Tutorial). Detailed documentation can be found [here](https://github.com/broadinstitute/PANOPLY/wiki).

### Contact

Email proteogenomics@broadinstitute.org with questions, comments or feedback.
