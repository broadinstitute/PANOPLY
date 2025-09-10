# ```panoply_select_all_pairs```

## Description
Internal module for collapsing a paired array with (potentially) missing values into a fully-defined array. Used internally by [panoply_unified_workflow](./Pipelines%3A-panoply_unified_workflow) to exclude missing omic datasets from analyses.


## Input

* ```pairs_input```: (`Array[Pair[String?,File?]]+`) paired array which may contain missing values


## Output

* ```pair_string```: (`Array[String]`) array containing the left-pair values for all complete pairs
* ```pair_file```: (`Array[File]`) array containing the right-pair values for all complete pairs
* ```pairs```: (`Array[Pair[String,File]]`) paired array containing values for all complete pairs
