# FragPipe DDA search and spectral library creation
**Version**: FragPipe 19.0, MSFragger 3.6, Philosopher 4.7.0, IonQuant 1.8.9

## `panoply_fragpipe` — workflow to run FragPipe pipeline (equivalent to GUI version)
### Inputs
- `fragpipe_workflow` (File): Google Bucket path to FragPipe workflow configuration (.workflow). This file can be saved from FragPipe GUI upon selecting the suitable parameters
- `database` (File): Google Bucket path to background proteome file (.fasta)
- `files_folder` (Directory): Google Bucket path to folder containing raw DDA files to be searched (ex: "gs://fc-3f59fceb-8ce2-4855-9198-d6f6527cd8af/Experiment_1/")

### Outputs
- `fragpipe_output` (File): zip containing all FragPipe and additional modules' outputs
- `fragpipe_processed_data` (File): zip containing intermediary files (e.g., .mzBIN) created by FragPipe

### Configuring a custom manifest file (e.g., grouping by experiments or replicates):
Paths to the raw files must start with `/root/fragpipe_fragpipe/data/` (this is how it's implemented in the Docker that the workflow runs on). It is recommended that you create this file with FragPipe GUI and edit the path prefix to fit this requirement.
```
/root/fragpipe/data/file1.mzML	exp 1	DDA
/root/fragpipe/data/file2.mzML	exp	2	DDA
/root/fragpipe/data/file3.mzML	exp	3	DDA
```
Once you save this file and upload to Google Bucket, you'd then configure the following parameter:
- `fragpipe_manifest` (File): Google Bucket path to manifest file (.fp-manifest)