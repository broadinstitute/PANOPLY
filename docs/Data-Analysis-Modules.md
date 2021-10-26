# Description of PANOPLY Data Analysis Modules

The following links provide documentation on each PANOPLY Data Analysis Module:

* [panoply_association](./Data-Analysis-Modules%3A-panoply_association)
* [panoply_blacksheep](./Data-Analysis-Modules%3A-panoply_blacksheep)
* [panoply_cmap_analysis](./Data-Analysis-Modules%3A-panoply_cmap_analysis)
* [panoply_cna_correlation](./Data-Analysis-Modules%3A-panoply_cna_correlation)
* [panoply_cons_clust](./Data-Analysis-Modules%3A-panoply_cons_clust)
* [panoply_immune_analysis](./Data-Analysis-Modules%3A-panoply_immune_analysis)
* [panoply_mo_nmf](./Data-Analysis-Modules%3A-panoply_mo_nmf)
* [panoply_rna_protein_correlation](./Data-Analysis-Modules%3A-panoply_rna_protein_correlation)
* [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea)
* [panoply_quilts](./Data-Analysis-Modules%3A-panoply_quilts)

Each data analysis module is run in its own (appropriately named) subdirectory that contains results and intermediate files generate during execution, along with code.

The following runtime parameters apply to each task, they can be changed but have the following defaults unless otherwise specified:

* ```memory```: (Int, default = 2) number of GB of memory
* ```disk_space```: (Int, default = 10) local disk SSD
* ```num_threads```: (Int, default = 1) number of cpu threads
* ```num_preemptions```: (Int, default = 0) number of preemptible attempts