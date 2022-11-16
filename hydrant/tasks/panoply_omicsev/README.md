## PANOPLY Implementation of OmicsEV

Find full OmicsEV documentation [here](https://github.com/bzhanglab/OmicsEV).

### Mandatory Inputs

**General mandatory inputs:**

* `STANDALONE` (String): Either `"true"` or `"false"`. Determines whether OmicsEV should run as part of `panoply_main` pipeline with inputs from `panoply_harmonize`, or if OmicsEV should run independently with .gct inputs. Preference is to run with `STANDALONE == "false"` if possible.
* `yaml_file` (File): Parameters yaml file from PANOPLY setup. Must contain gene.id.col and groups information.
* `label` (String): Label for this data set (e.g. "BRCA-proteome"). HTML report will include this label in the file name. Tarball output will also be in this folder if `STANDALONE == "false"`.

**Mandatory inputs if `STANDALONE = "false"`:**

* `panoply_harmonize_tar` (File): Tar output from panoply_harmonize module.

**Mandatory inputs if `STANDALONE = "true"`:**

* `data_files` (Array[File]): Array of .gct files to be used for OmicsEV. These are typically all the same -ome from the same samples with some variation in how abundance values were calculated. For example, multiple proteome files generated using different batch correction methods. Array[File] format is `["/path/to/file-1.gct","/path/to/file-2.gct","/path/to/file-3.gct"]`. *Make sure there is no space between each file, just a single comma*. A single data file input is acceptable (e.g. `["/path/to/file-1.gct"]`).
* `sample_anno_file` (File): A .csv file with sample annotations. File must include same samples as are in `data_files`.

### Optional Inputs

**General optional inputs:**

* `class_column_name` (String): Column in sample annotation file that contains class information. This should be a categorical annotation (e.g. sex, cancer type, etc.) to assess differences between phenotypes.
* `batch_column_name` (String): Column in sample annotation file that contains batch information. This is used to assess batch effects. Must be integer values with no NA's. If `STANDALONE = true`, this can also be a column in the column-description of each .gct file in `data_files`.
* `data_log_transformed` (Boolean): Indicates whether the protein-level data has been log transformed.
* `rna_log_transformed` (Boolean): Indicates whether the rna data has been log transformed.
* `do_function_prediction` (Boolean): Whether or not to perform gene function prediction analysis. If this is not necessary, it is suggested to turn this off. Gene function prediction takes a considerable amount of time. Default is `true`.

**Optional inputs if `STANDALONE = true`:**

* `rna_file` (File): A .gct file of RNAseq data. Used only for RNA-protein correlation at the end of OmicsEV.
* `data_type` (String): Either "protein" or "gene". Describes the data type of `data_files`. Default is "protein".

### Output

* `report` (File): An HTML report of OmicsEV evaluation metrics for each dataset.
* `outputs` (File): A tarball with the full OmicsEV outputs, including intermediary tables and figures. 
