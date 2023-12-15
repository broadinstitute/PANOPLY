# msConvert DIA RAW file Conversion to MZML
**Version**: msConvert 3.0

## Description
This workflow wraps the msConvert GUI and converts DIA RAW files into MZML outputs. MZML outputs are moved to a specified google bucket to allow for subsequent analysis.

## `panoply_msconvert` — workflow to run msConvert
### Inputs
**Workflow inputs**
- `files_folder` (Directory): Google Bucket path to folder containing raw files to be searched (ex: "gs://fc-3f59fceb-8ce2-4855-9198-d6f6527cd8af/Experiment_1/")
- `output_folder` (String): Google Bucket path to folder where MZML files will localize folloing completion of the workflow
- `file_of_files` (File, optional): Google Bucket path to a text file storing Google Bucket paths to raw files to be analyzed (only works with Thermo data as timsTOF `.d` are technically folders); pecifying this would overwrite `files_folder`
- `staggered_window` (Boolean, defalut="False"): for DIA RAW files with 50% overlapping isolation windows
- `massError` (String, default="10.0"): demultiplex parameter for staggard window conversion

**Terra parameters**
- `local_disk_gb` (Integer, default=1000): GB of storage space in the compute instance
- `num_cpus` (Integer, default=32): number of cores in the compute instance
- `num_preemptions` (Integer, default=0): number 
- `ram_gb` (Integer, default=128): GB of RAM memory in the compute instance

## Developer guide 
### Flow
1. Dockers are sourced from third parties and should remain current with updates
2. Files are localized into `/root/raw_files/data` data folder from `files_folder` or `file_of_files` (if supplied)
3. (Experimental) `try_one_file` flag removes all but one file in data folder. Used to test if workflow works without running on the whole dataset

### Paths
- Working directory is `/root/raw_files`
- Input data will be in `/root/raw_files/data`
- Outputs will be in user specified output directory

## References
- Chambers, M. C., Maclean, B., Burke, R., Amodei, D., Ruderman, D. L., Neumann, S., Gatto, L., Fischer, B., Pratt, B., Egertson, J., Hoff, K., Kessner, D., Tasman, N., Shulman, N., Frewen, B., Baker, T. A., Brusniak, M.-Y., Paulse, C., Creasy, D., … Mallick, P. (2012). A cross-platform toolkit for mass spectrometry and proteomics. Nature Biotechnology, 30(10), Article 10.

