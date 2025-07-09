# Spectronaut SNE Merge and Combine
**Version**: Spectronaut v20, cloud enabled

## Description
This workflow invokes the Spectronaut CLI interface. The workflow is made for SNE combination and management using the merge function. More detailed documentation can be found in the [Spectronaut 20 User Manual](https://github.com/broadinstitute/PANOPLY/blob/dev/third-party-modules/panoply_spectronaut/Spectronaut-20-Manual.pdf)

## `panoply_spectronaut` — workflow to run Spectronaut on the cloud
### Inputs
**Workflow inputs**
- `experiment_name` (String): Name of the experiment used by Spectronaut. Outout files are written in a directory with `experiment_name`.
- `fasta` (File): Google Bucket path to protein search database (.fasta) (ex: "gs://fc-7579f17f-822d-47be-af67-5c99fc597005/20220222uniprotproteomeUP000005640.bgsfasta")
- `files_folder` (Directory): Google Bucket path to folder containing SNE files to be merged or combined (ex: "gs://fc-7579f17f-822d-47be-af67-5c99fc597005/RawFiles/Ultra")
- `sne_out`  (Boolean): Parameter determining if output is an SNE file, default is False using the combine funciton, True uses merge function
- `settings_schema` (File, optional): Google Bucket path to a settings schema file exported from the Spectronaut UI, used for combine function
- `fasta_1` (File): Google Bucket path to additional protein search database (.fasta)
- `report_schema` (File, optional): Spectronaut report scheme for specifying specific output report formats
- `report_schema_1` (File, optional): Additional Spectronaut report scheme for specifying specific output report formats
- `report_schema_2` (File, optional): Additional Spectronaut report scheme for specifying specific output report formats
- `file_of_files` (File, optional): Google Bucket path to a text file storing Google Bucket paths to raw files to be analyzed (only works with Thermo data as timsTOF `.d` 

**Terra parameters**
- `local_disk_gb` (Integer, default=1000): GB of storage space in the compute instance
- `num_cpus` (Integer, default=32): number of cores in the compute instance
- `num_preemptions` (Integer, default=0): number 
- `ram_gb` (Integer, default=128): GB of RAM memory in the compute instance

### Outputs
- `spectronaut_output.zip` (File):
    - All Spectronaut outputs are contained in the output zip file.

