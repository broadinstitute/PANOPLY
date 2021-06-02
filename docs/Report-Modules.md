# Description of PANOPLY Report Modules

PANOPLY Report Modules generate an interactive R markdown report rendered into HTML to summarize the results generated from the corresponding Analysis Module. The following links provide documentation on each PANOPLY report module:

* [panoply_association_report](./Report-Modules%3A-panoply_association_report)
* [panoply_blacksheep_report](./Report-Modules%3A-panoply_blacksheep_report)
* [panoply_cna_correlation_report](./Report-Modules%3A-panoply_cna_correlation_report)
* [panoply_cons_clust_report](./Report-Modules%3A-panoply_cons_clust_report)
* [panoply_immune_analysis_report](./Report-Modules%3A-panoply_immune_analysis_report)
* [panoply_mo_nmf_report](./Report-Modules%3A-panoply_mo_nmf_report)
* [panoply_normalize_ms_data_report](./Report-Modules%3A-panoply_normalize_ms_data_report)
* [panoply_rna_protein_correlation_report](./Report-Modules%3A-panoply_rna_protein_correlation_report)
* [panoply_sampleqc_report](./Report-Modules%3A-panoply_sampleqc_report)
* [panoply_ssgsea_report](./Report-Modules%3A-panoply_ssgsea_report)

The following runtime parameters apply to each task, they can be changed but have the following defaults unless otherwise specified:

* ```memory```: (Int, default = 2) number of GB of memory
* ```disk_space```: (Int, default = 10) local disk SSD
* ```num_threads```: (Int, default = 1) number of cpu threads
* ```num_preemptions```: (Int, default = 0) number of preemptible attempts