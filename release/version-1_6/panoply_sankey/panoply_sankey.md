# ```panoply_sankey```

## Description

PANOPLY workflow to compare sample-membership within a given annotation, between 2+ annotation-files. For each pair of annotation-files, a sankey diagram will be created comparing the sample-membership in the first file to those in the second file. If a "primary" annotation file is provided, the module will additionally create 3-way sankey diagrams, with the values from the primary-annotation-value in the center.

## Input

### Required inputs:
* ```annot_files```: (Array[File]+) array of sample-annotation files to be compared (accepts multiple file-format, but prefers `.tsv`)
* ```annot_file_labels```: (Array[String]+) array of labels uniquely identifying each annotation file (must match length and order of `annot_files`)
* ```label```: (String) prefix for naming the output tar file

### Optional inputs:
* ```annot_file_primary```: (File) sample-annotation file which should be highlighted/centered in diagrams (accepts multiple file-format, but prefers `.tsv`)
* ```annot_label_primary```: (String) label uniquely identifying the primary annotation file

* ```id_column```: id column for identifying entries (e.g. "Sample.ID"); uses rownames if not provided
* ```annot_column```: annotation column for sankey comparisons (e.g. "NMF.consensus")
* ```annot_prefix```: prefix to prepend to annotation values (e.g. "C" -> C1 C2 C3, instead of 1 2 3)
* ```color_str```: (String) comma-separated string of colors to use for sankey diagrams (e.g. '#AA0000,#0000AA' will assign annotations-subvalues colors on a gradient from red to blue)
 

## Output

* ```tar_out```: (`.tar` file) tarball containing HTML and PDF copies of all sankey diagrams