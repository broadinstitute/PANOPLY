# Spectronaut search and spectral library creation
**Version**: Spectronaut v18, cloud enabled

## Description
This workflow invokes the Spectronaut CLI interface. All of the settings can be configured in Spectronaut UI and saved as a `.prop` file to be uploaded to this workflow. The settings can also be edited using a text editor. More detailed documentation can be found in the Spectronaut 18 User Manual

## `panoply_spectronaut` — workflow to run Spectronaut on the cloud
### Inputs
**Workflow inputs**
- `experiment_name` (String): Name of the experiment used by Spectronaut. Outout files are written in a directory with `experiment_name`.
- `fasta` (File): Google Bucket path to protein search database (.fasta) (ex: "gs://fc-7579f17f-822d-47be-af67-5c99fc597005/20220222uniprotproteomeUP000005640.bgsfasta")
- `files_folder` (Directory): Google Bucket path to folder containing raw data files to be searched (ex: "gs://fc-7579f17f-822d-47be-af67-5c99fc597005/RawFiles/Ultra")
- `analysis_settings` (File,optional): Google Bucket path to Spectronaut settings file (.prop). This file can be saved from Spectronaut UI upon selecting the suitable parameters, or default settings can be obtained from a Spectronaut/Biognosys
- `condition_setup` (File, optional): Google Bucket path to a `.tsv` file exported from the Spectronaut UI
- `file_of_files` (File, optional): Google Bucket path to a text file storing Google Bucket paths to raw files to be analyzed (only works with Thermo data as timsTOF `.d` are technically folders); specifying this would overwrite `files_folder`
- `json_settings` (File,optional): Google Bucket path to a `json` file with Spectronaut settings. This is an alternate way of specifying settings (instead of a `.prop` file)
- `fasta_1` (File): Google Bucket path to additional protein search database (.fasta)
- `license_key` (String, optional): The license key needed to activate and run Spectronaut. A default key is automatically provided. This input can be used to override the default key, or if the default key is invalid or expired
- `report_scheme` (File, optional): Spectronaut report scheme for specifying specific output report formats
- `spectral_library` (File, optional): Spectral library to be used in DIA search. If a spectral library is not specified, Spectronaut runs a Direct DIA search
- `spectral_library_1` (File, optional): Additional spectral library to be used in DIA search
- `pulsar_settings` (File, optional): Google Bucket path to a Pulsar settings file. When running Spectronaut spectral library generation, this file is used to configure Pulsar 

**Terra parameters**
- `local_disk_gb` (Integer, default=1000): GB of storage space in the compute instance
- `num_cpus` (Integer, default=32): number of cores in the compute instance
- `num_preemptions` (Integer, default=0): number 
- `ram_gb` (Integer, default=128): GB of RAM memory in the compute instance

### Outputs
- `spectronaut_output.zip` (File):
    - All Spectronaut outputs are contained in the output zip file.

