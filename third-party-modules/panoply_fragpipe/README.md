# FragPipe DDA search and spectral library creation
**Version**: FragPipe 19.1, MSFragger 3.7, Philosopher 4.8.0, IonQuant 1.8.10

## Description
This workflow wraps the FragPipe program using the CLI interface. All of the configurations can be configured in FragPipe GUI and saved as a `.workflow` file to be uploaded to this workflow. It is also capable of generating a spectral library via EasyPQP. Outputs are separated into converted files and generated files with final ion, PSM, protein level report tables.

Step-by-step slides for running this workflow are available [here](https://docs.google.com/presentation/d/1Rmer-oOaP-Eqo5tyxtGG-AlM4BpQaWjkpL0GBxd-VlY/edit?usp=sharing).

## `panoply_fragpipe` — workflow to run FragPipe pipeline (equivalent to GUI version)
### Inputs
**Workflow inputs**
- `database` (File): Google Bucket path to background proteome file (.fasta)
- `files_folder` (Directory): Google Bucket path to folder containing raw DDA files to be searched (ex: "gs://fc-3f59fceb-8ce2-4855-9198-d6f6527cd8af/Experiment_1/")
- `fragpipe_workflow` (File): Google Bucket path to FragPipe workflow configuration (.workflow). This file can be saved from FragPipe GUI upon selecting the suitable parameters
- `file_of_files` (File, optional): Google Bucket path to a text file storing Google Bucket paths to raw files to be analyzed (only works with Thermo data as timsTOF `.d` are technically folders); pecifying this would overwrite `files_folder`
- `fragpipe_manifest` (File): Google Bucket path to FragPipe manifest file (.fp-manifest)
- `raw_file_type` (String, default="DDA"): acquisition type of raw data, needed for automatic creation of the manifest

**Terra parameters**
- `local_disk_gb` (Integer, default=1000): GB of storage space in the compute instance
- `num_cpus` (Integer, default=32): number of cores in the compute instance
- `num_preemptions` (Integer, default=0): number 
- `ram_gb` (Integer, default=128): GB of RAM memory in the compute instance

### Outputs
- `fragpipe_output.zip` (File):
    - All FragPipe long-format reports of properties, metrics, and intensities (`library.tsv`, `psm.tsv`, `ion.tsv`, `protein.tsv`)
    - Search parameter files (`fragger.params`, `msbooster_params.txt`, ...)
    - Intermediary files (`.pepXML`, `.pin`, `.peakpkl`, ...)
    - Outputs of other modules if specified to run (e.g., `.speclib` from EasyPQP)
    - Log files for main pipeline and specific modules (`log.txt`, `filter.log`)
- `fragpipe_processed_data.zip` (File):
    - Input files used in the search
    - Converted files (if created, `.mzML`) created by FragPipe

### Configuring a custom manifest file (e.g., grouping by experiments or replicates):
Paths to the raw files must start with `/root/fragpipe_fragpipe/data/` (this is how it's implemented in the Docker that the workflow runs on). It is recommended that you create this file with FragPipe GUI and edit the path prefix to fit this requirement.
```
/root/fragpipe/data/file1.mzML  exp 1   DDA
/root/fragpipe/data/file2.mzML  exp 2   DDA
/root/fragpipe/data/file3.mzML  exp 3   DDA
```
Once you save this file and upload to Google Bucket, you'd then configure the following parameter:
- `fragpipe_manifest` (File): Google Bucket path to manifest file (.fp-manifest)

## Developer guide 
### Flow
1. Docker contains MSFragger, philosopher, IonQuant, easyPQP that FragPipe headless uses
2. Files are localized into `/root/fragpipe/data` data folder from `files_folder` or `file_of_files` (if supplied)
3. Inputted database will be written into the workflow file before running FragPipe
4. (Experimental) `try_one_file` flag removes all but one file in data folder. Used to test if workflow works without running on the whole dataset
5. Manifest can either be user-supplied or automatically generated with `get_fp_manifest.py` script which doesn't assume any experimental design
6. Call FragPipe CLI command with supplied workflow, localized files, and manifest (auto-generated or supplied)

### Paths
- Working directory is `/root/fragpipe`
- Input data will be in `/root/fragpipe/data`
- Outputs will be in `/root/fragpipe/out`

## References
- Kong, A. T., Leprevost, F. V., Avtonomov, D. M., Mellacheruvu, D., & Nesvizhskii, A. I. (2017). MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry–based proteomics. Nature Methods, 14(5), 513-520.
- Yu, F., Teo, G. C., Kong, A. T., Haynes, S. E., Avtonomov, D. M., Geiszler, D. J., & Nesvizhskii, A. I. (2020). Identification of modified peptides using localization-aware open search. Nature Communications, 11(1), 1-9.
- Polasky, D. A., Yu, F., Teo, G. C., & Nesvizhskii, A. I. (2020). Fast and Comprehensive N-and O-glycoproteomics analysis with MSFragger-Glyco. Nature Methods, 17, 1125-1132.
- Chang, H. Y., Kong, A. T., da Veiga Leprevost, F., Avtonomov, D. M., Haynes, S. E., & Nesvizhskii, A. I. (2020). Crystal-C: A computational tool for refinement of open search results. Journal of Proteome Research, 19(6), 2511-2515.
- Geiszler, D. J., Kong, A. T., Avtonomov, D. M., Yu, F., da Veiga Leprevost, F., & Nesvizhskii, A. I. (2020). PTM-Shepherd: analysis and summarization of post-translational and chemical modifications from open search results. Molecular & Cellular Proteomics.
- da Veiga Leprevost, F., Haynes, S. E., Avtonomov, D. M., Chang, H. Y., Shanmugam, A. K., Mellacheruvu, D., Kong, A. T., & Nesvizhskii, A. I. (2020). Philosopher: a versatile toolkit for shotgun proteomics data analysis. Nature Methods, 17(9), 869-870.
- Yu, F., Haynes, S. E., Teo, G. C., Avtonomov, D. M., Polasky, D. A., & Nesvizhskii, A. I. (2020). Fast quantitative analysis of timsTOF PASEF data with MSFragger and IonQuant. Molecular & Cellular Proteomics.
- Yu, F., Haynes, S. E., & Nesvizhskii, A. I. (2021). IonQuant enables accurate and sensitive label-free quantification with FDR-controlled match-between-runs. Molecular & Cellular Proteomics, 20.
- Teo, G. C., Polasky, D. A., Yu, F., Nesvizhskii, A. I. (2020). A fast deisotoping algorithm and its implementation in the MSFragger search engine. Journal of Proteome Research.
- Tsou, C. C., Avtonomov, D., Larsen, B., Tucholska, M., Choi, H., Gingras, A. - C., & Nesvizhskii, A. I. (2015). DIA-Umpire: comprehensive computational framework for data-independent acquisition proteomics. Nature methods, 12(3), 258-264.
