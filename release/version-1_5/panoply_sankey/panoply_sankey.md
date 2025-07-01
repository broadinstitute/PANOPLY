# ```panoply_sankey_workflow```

## Description

PANOPLY workflow to compare sample-membership in a given annotation, between two or more annotation-files. Can optionally take a "primary" annotation file, to be highlighted in the report, and placed in the center of all three-way sankey-diagrams.  It is used in the [panoply_nmf_workflow](https://github.com/broadinstitute/PANOPLY/wiki/Pipelines%3A-panoply_nmf_workflow) workflow compare the results of unsupervised non-negative matrix factorization (NMF)-based clustering, across different clustering-attempts. 

The workflow expects an array of sample-annotation file with shared sample-IDs, and a column specifying which annotation-values should be visualized.


The workflow consists of two modules: 

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_sankey</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_sankey)         |   create sankey diagram comparison figures |
| [<code>panoply_sankey_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_sankey_report)          |   RMarkown report for `panoply_sankey`  |



## Input

### Required inputs:
* ```annot_files```: (Array[File]+) array of sample-annotation files to be compared (accepts multiple file-format, but prefers `.tsv`)
* ```annot_labels```: (Array[String]+) array of labels uniquely identifying each annotation file (must match length and order of `annot_files`)


### Optional inputs:
* ```annot_file_primary```: (File) sample-annotation file which should be highlighted/centered in diagrams (accepts multiple file-format, but prefers `.tsv`)
* ```annot_label_primary```: (String) label uniquely identifying the primary annotation file

* ```id_column```: id column for identifying entries (e.g. "Sample.ID"); uses rownames if not provided
* ```annot_column```: annotation column for sankey comparisons (e.g. "NMF.consensus")
* ```annot_prefix```: prefix to prepend to annotation values (e.g. "C" -> C1 C2 C3, instead of 1 2 3)
 
 
## Output

* ```sankey tar - ${label}_sankey_diagrams.tar```: Sankey diagrams comparisons between the NMF-clustering-results of supplied -omes.
* ```sankey report - ${label}_sankey_rmd.html```: Summary report of sankey diagram comparisons between the NMF-clustering-results of supplied -omes.

