The compare.R script can be used to compare the results of two PANOPLY runs against each other.

Usage:
  Rscript compare.R --run1 path/to/run1 --run2 path/to/run2 [-o report.txt]


The script relies on three helper files:
  * `helpers.R` - This script contains various helper functions, used to normalize filenames and compare file contents across runs. It provides support for comparing files in a number of ways (e.g. strictly, within a threshold, across a subset of rows, etc.)
  * `module_checks.R` - This script contains custom comparison functions, for use with certain stochastic modules that require a more nuanced comparison method (e.g. panoply_nmf).
  * `task_rules.yaml` - This config file provides a ruleset for each module, describing what ruleset should be used to compare output files. When a new module is added to PANOPLY, a new ruleset for that module should be added to this file.
  * `folder_map.yaml` - This config file describes the file-structure of the unified pipeline assembly; it provides a ruleset for how files in a directory should be compared.

This folder also contains a development tool, `validate_rules.R`, which can be used to validate the `rules.yaml` file against the current PANOPLY pipeline. The validation script checks for modules which are missing from the current rules file and need to be added.