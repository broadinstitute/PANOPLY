## Introduction

Welcome to the PANOPLY wiki!

[PANOPLY](https://github.com/broadinstitute/PANOPLY) is a cloud-based platform for applying statistical and machine learning algorithms to transform multi-omic data from cancer samples into biologically meaningful and interpretable results. PANOPLY leverages [Terra](http://app.terra.bio), and is designed to be flexible, automated, reproducible, scalable, and secure. A wide array of [algorithms](https://github.com/broadinstitute/PANOPLY/wiki) for proteogenomic data analysis have been implemented in PANOPLY.

Documentation and Tutorial can be found in the **Sidebar**.

## Requirements

1. To run PANOPLY, you need a Terra account. To create an account follow the steps [here](https://support.terra.bio/hc/en-us/articles/360034677651-Account-setup-and-exploring-Terra). A Google account is required to create a Terra account.
2. Once you have an account, you need to create a [cloud billing account and project](https://support.terra.bio/hc/en-us/articles/360026182251-How-to-set-up-billing-in-Terra). Here, you will find step-by-step instructions to set up billing, including a free credits program for new users to try out Terra with $300 in Google Cloud credits.

## Inputs

Input data consists of a **ZIP** compressed file containing the following data tables:

##### Required:
* Datasets for [*preconfigured piplelines*](https://app.terra.bio/#workspaces/broad-firecloud-cptac/PANOPLY_Production_Pipelines_v1_0):
    - At least **one proteomics dataset** (global proteome, phosphoproteome, acetylome, ubiquitylome) --- after database searching and quantification;
    - Genomics data --- CNA and RNA data (both are required), normalized and filtered as needed;
* Dataset(s) for running [*individual analysis modules*](https://app.terra.bio/#workspaces/broad-firecloud-cptac/PANOPLY_Production_Modules_v1_0):
    - At least **one dataset**: either proteomics (global proteome, phosphoproteome, acetylome, ubiquitylome--after database searching and quantification) or genomics (RNA, CNA--normalized and filtered as needed) data;
* `Annotation` table with at least `Sample.ID` and `Type` columns, in addition to an arbitrary number of annotations for the samples.

##### Optional:
* A `groups` file that lists a subset of annotations (one per line) to be used for association and enrichment analysis;
* A default `parameter` file is used, unless one is specified in the input (in `yaml` format);
* Pathway databases for PTM-SEA (PTM-SigDB v1.9.0) and GSEA (v6.2 hallmark pathways) are automatically included, and can be over-ridden by providing appropriate `gmt` input files.

##### Input File Specifications:
|File    |Format         | Specification                                                                                |
|--------|---------------|----------------------------------------------------------------------------------------------|
|`annotation`    | `csv` | The `annotation` table must include `Sample.ID` and `Type` columns. `Sample.ID`s must be unique and cannot have duplicates. If replicate samples are present, they must have unique `Sample.ID`s, but can be identified as replicates by using an additional `Participant` annotation column with identical id's for replicate samples. The `annotation` table input to the `panoply_parse_sm_table` module must include `Experiment` (TMT or iTRAQ plex number) and `Channel` (TMT or iTRAQ channel id). |
|`groups`        | `csv` | Lists annotation columns (one per line) present in the `annotation` table for use in association and enrichment analysis in various PANOPLY modules. If a `groups` file is not specified in the input, one will be interactively created by the `PANOPLY-startup-notebook` using the annotations present in the `annotation` table. |
|proteomics data | [`gct` v1.3](https://clue.io/connectopedia/gct_format) | Can be proteome, phosphoproteome, acetylome, ubiquitylome. Data matrix (resulting from MS database search) must have protein/PTM-site as rows and samples as columns, with entries containing abundance values (usually log2 ratio to a common reference). At least **one** proteomics data type must be specified. PANOPLY provides an option to normalize proteomics data, if needed. |
|genomics data   | [`gct` v1.3](https://clue.io/connectopedia/gct_format) | **Both** CNA (normalized log-ratio, derived from WXS, WGS or combination) and RNA expression (log-transformed and normalized, derived from RNAseq) data are required. These data *must* be normalized prior to input in PANOPLY. To include *mutation data*, samples with mutations in relevant genes must be identified and an annotation column added for each gene in the `annotation` table---`maf` or `vcf` cannot be directly input. |
|parameters     | [`yaml`](https://yaml.org/) | The [`master-parameters.yaml` file](https://github.com/broadinstitute/PANOPLY/blob/release-1_0/panda/panda-src/defaults/master-parameters.yaml) (also available in the PANOPLY Production Workspaces on Terra) contains default parameter settings for all PANOPLY modules, and is used as a default, if no input parameters file is specified. |
|pathway databases | [`gmt`](https://www.genepattern.org/file-formats-guide#GMT) | If no pathway databases are specified in the input, [PTM-SigDB v1.9.0](https://github.com/broadinstitute/PANOPLY/blob/release-1_0/panda/panda-src/defaults/ptm.sig.db.all.uniprot.human.v1.9.0.gmt) is used for PTM-SEA, and the [Cancer Hallmark Pathways v6.2](https://github.com/broadinstitute/PANOPLY/blob/release-1_0/panda/panda-src/defaults/h.all.v6.2.symbols.gmt) is used for GSEA and ssGSEA.

##### Notes:

* All proteomics and genomics data tables must have sample ids (column names) conforming to the `Sample.ID`s used in the `annotation` table. Ideally, the sample annotations included in the `annotation` table would be present in the GCT v1.3 data files as column annotations.
* All genomics data tables must be appropriately normalized/filtered prior to use in PANOPLY. Proteomics data can be optionally normalized/filtered in PANOPLY.


## Documentation

PANOPLY modules are grouped into `Data Preparation`, `Data Analysis`, `Report` and `Support` categories. There is also documentation on navigating results. In addition, a step-by-step tutorial is included. The documentation main pages are listed below, and can also be accessed using the sidebar for this Wiki.

* [Tutorial](https://github.com/broadinstitute/PANOPLY/wiki/PANOPLY-Tutorial)
* [Pipelines](https://github.com/broadinstitute/PANOPLY/wiki/Pipelines)
* [Data Preparation Modules](https://github.com/broadinstitute/PANOPLY/wiki/Data-Preparation-Modules)
* [Data Analysis Modules](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules)
* [Report Modules](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules)
* [Support Modules](https://github.com/broadinstitute/PANOPLY/wiki/Support-Modules)
* [Navigating Results](https://github.com/broadinstitute/PANOPLY/wiki/Navigating-Results)
* [PANOPLY without Terra](https://github.com/broadinstitute/PANOPLY/wiki/PANOPLY-without-Terra)
* [Customizing PANOPLY](https://github.com/broadinstitute/PANOPLY/wiki/Customizing-PANOPLY)

