# Description of PANOPLY Data Preparation Modules

The following links provide documentation on each PANOPLY Data Preparation Module:

* [panoply_normalize_ms_data](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)
* [panoply_parse_sm_table](./Data-Preparation-Modules%3A-panoply_parse_sm_table)
* [panoply_sampleqc](./Data-Preparation-Modules%3A-panoply_sampleqc)

Each data preparation module is run in its own (appropriately named) subdirectory that contains results and intermediate files generate during execution, along with code.

The following runtime parameters apply to each task, they can be changed but have the following defaults unless otherwise specified:

* ```memory```: (Int, default = 2) number of GB of memory
* ```disk_space```: (Int, default = 10) local disk SSD
* ```num_threads```: (Int, default = 1) number of cpu threads
* ```num_preemptions```: (Int, default = 0) number of preemptible attempts
