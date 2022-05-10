# ```panoply_harmonize```

## Description

This module harmonizes RNA, CNA, and proteome data to create gene-centric tables with common samples as columns and common genes as rows.

## Input

Required inputs:

* ```inputData```: (`.tar` file) tarball from ```panoply_rna_protein_correlation```
* ```rnaExpr```: (`.gct` file) RNA expression data
* ```cnaExpr```: (`.gct` file) CNA expression data
* ```type```: (String) data type
* ```standalone```: (String) 
* ```yaml```: (`.yaml` file) master-parameters.yaml

Optional inputs:

* ```analysisDir```: (String) name of analysis directory
* ```pomeGeneIdCol```: (String, default = 'GeneSymbol') gene id column in proteome data
* ```cnaGeneIdCol```: (String, default = 'id') gene id column in CNA expression data
* ```rnaGeneIdCol```: (String, default = 'id') gene id column in RNA expression data
* ```ndigits```: (Int, default = 5) number of decimal digits to use in output tables
* ```na_max```: (Float, default = 0.7) maximum allowed NA values per protein/site/row; can be fraction.
* ```duplicate_gene_policy```: (String, default = 'maxvar')
* ```gene_id_col```: (String, default = 'geneSymbol') name of sample annotation column containing gene ids.
* ```outFile```: (String, default = "panoply_harmonize-output.tar") output `.tar` file name


## Output

Tarball including the following files in the `harmonized-data` subdirectory:

* Harmonized data matrices for CNA (`cna-matrix.csv`), RNA (`rna-matrix.csv`) and proteome (`*-matrix.csv`) with samples (columns) and genes (rows) in identical order in all tables.
* Matching sample annotation for samples included in the data matrices (`sample-info.csv`)
* Class vectors files (`*.cls`) for all sample annotations
