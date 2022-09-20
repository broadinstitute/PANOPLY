## PANOPLY implementation of COSMO

Find full COSMO documentation [here](https://github.com/bzhanglab/COSMO).

### Inputs

Mandatory inputs:

* `STANDALONE` (Boolean): Determines whether COSMO should run as part of `panoply_main` pipeline with inputs from `panoply_harmonize`, or if COSMO should run independently with .gct inputs.
* `yaml_file` (File): Parameters yaml file from PANOPLY setup. Must contain gene.id.col and groups information.

Mandatory inputs if `STANDALONE = false`:

* `panoply_harmonize_tar` (File): Tar output from panoply_harmonize module.

Mandatory inputs if `STANDALONE = true`:

* `d1_file` (File): A .gct file. Typically proteome or other protein-level data.
* `d2_file` (File): A .gct file. Typically RNAseq file to be compared with `d1_file` for sample mislabeling.
* `sample_file` (File): A .csv file with sample annotations.

Optional inputs:

* `sample_label` (String): Column header(s) in `sample_file` to use for prediction. This typically includes 'gender'. Column will be excluded if it is unbalanced (likely to cause an error). If multiple inputs, separate inputs with a comma (e.g. 'gender,msi'). If no input, columns are pulled from groups in yaml file.

### Output

The main output is `output.zip`, a zip file containing the key outputs from the COSMO pipeline. This can be fed into `panoply_cosmo_report` for report generation.

### Common Pitfalls

The main source of errors in COSMO is improper selection of sample labels. PANOPLY-specific preprocessing should eliminate sample labels that are likely to cause errors. COSMO requires clinical attributes that are well-balanced with only two levels (e.g. male/female, positive/negative). These are used to predict if there is mislabeling between the sample annotation file (`sample_file`) and any of the data files (`d1_file`, `d2_file`). 

If you get the following error message, this likely means that one of the sample label columns in not well-balanced. Consider removing it as an input. 

```<simpleError in { which = foldid == i if (length(dim(y)) > 1) y_sub = y[!which, ] else y_sub = y[!which] if (is.offset) offset_sub = as.matrix(offset)[!which, ] else offset_sub = NULL glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda, offset = offset_sub, weights = weights[!which], ...)}: task 4 failed - "one multinomial or binomial class has 1 or 0 observations; not allowed">```
