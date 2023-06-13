
# ```panoply_so_nmf_assemble_results```

## Description

This module assembles the outputs of multiple [panoply_so_nmf_gct](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_so_nmf) module-runs into a single results tar and a single reports tar.

## Input

### Required inputs:

* ```so_nmf_tar```: (Array[File]) Tar file output(s) of [<code>panoply_mo_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf) applied to single-ome(s)
* ```so_nmf_report```: (Array[File]) Tar report output(s) of [<code>panoply_mo_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_mo_nmf_report) applied to single-ome(s)
* ```so_nmf_ssgsea_tar```: (Array[File]) Tar file output(s) of [<code>panoply_ssgsea</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea) applied to single-ome(s)
* ```so_nmf_ssgsea_report```: (Array[File]) Tar file output(s) of [<code>panoply_ssgsea_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_ssgsea_report) applied to single-ome(s)


## Output

* ```so_nmf_resuls.tar```: Tarred results of NMF clustering and SSGSEA analysis, for each of the supplied -omes.
* ```so_nmf_reports.tar```: Tarred reports for NMF clustering and SSGSEA analysis, for each of the supplied -omes.
  