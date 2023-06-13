# ```panoply_so_nmf_sankey_workflow```

## Description

PANOPLY workflow to compare results of unsupervised non-negative matrix factorization (NMF)-based clustering, on a group of single -ome datasets. Creates sankey diagrams comparing the sample-composition of clusters, bewteen each single-ome NMF result. Can optionally take the results of multi-omics clustering, to compare the sample-composition of clusters using features from *all* -omes, to clusters using features from a single-ome.

The workflow expects the tarred output of [<code>panoply_so_nmf_assemble_results</code>](https://github.com/broadinstitute/PANOPLY/wiki/Support-Modules%3A-panoply_so_nmf_assemble_results). The workflow can optionally take the tarred output of [multi-omics NMF-based clustering](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf), for comparison between single-omic and multi-omic results.


The workflow consists of two modules: 

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_so_nmf_sankey</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_so_nmf_sankey)         |   create sankey diagram comparisons of so-NMF clustering-results |
| [<code>panoply_so_nmf_sankey_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_so_nmf_sankey_report)          |   RMarkown report for `panoply_so_nmf_sankey`  |



## Input

### Required inputs:
* ```so_nmf_tar```: (`.tar` file) Tar file output of [assembled single-ome NMF-based clustering results](https://github.com/broadinstitute/PANOPLY/wiki/Support-Modules%3A-panoply_so_nmf_assemble_results)
* ```label```: (String) label

### Optional inputs:
* ```mo_nmf_tar```: (`.tar` file)  Tar file output of [multi-omics NMF-based clustering](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf)
 
 
 
## Output

* ```sankey tar - ${label}_sankey_diagrams.tar```: Sankey diagrams comparisons between the NMF-clustering-results of supplied -omes.
* ```sankey report - ${label}_so_nmf_rmd.html```: Summary report of sankey diagram comparisons between the NMF-clustering-results of supplied -omes.

