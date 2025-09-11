# ```panoply_clumps_ptm_diffexp```

## Description

This module performs differential-expression analysis across all features from the provided PTM datasets, for each of the annotations provided in `groupsFile`. It produces `.tsv` file(s) for each annotation-column, which can be used as input to [Clumps-PTM](./Data-Analysis-Modules%3A-panoply_clumps_ptm) analysis.

## Input

### Required inputs:

* ```pSTY_gct```: (`.gct` file) phosphoproteome data matrix
* ```acK_gct```: (`.gct` file) acetylome data matrix
* ```ubK_gct```: (`.gct` file) ubiquitylome data matrix

* ```groupsFile```: (`.csv` file) annotation file, subsetted to annotations of interest for this analysis

* ```output_prefix```: (String, default="results") prefix used to name the output tar file
* ```yaml_file```: (`.yaml` file) master-parameters.yaml

### Optional inputs:

* ```sample_id_col```: (String, default="Sample.ID") column in `groupsFile` with sample identifiers.
* ```fdr_cutoff```: (Float, default=`NULL`) FDR threshold for significant features. Set to NULL to keep all features regardless of significance.
* ```min_samples```: (Int) minimum number of samples an annotation-subgroup can have before it is dropped from analysis (5 recommended, 3 absolute minimum)
* ```max_annot_levels```: (Int) maximum number of levels an annotation can have and be considered discrete

## Output

* ```diff_exp_files```: `.tsv` file(s) containing results of differential-expression analysis, for each of the annotation-columns in `groupsFile`.
* ```log_file```: (`.csv` file) Log file which lists every annotation tested, along with whether analysis was valid and had results. 