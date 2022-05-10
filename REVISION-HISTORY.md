### PANOPLY Revision History

**v1.1**

* New Modules:
  - ```panoply_mimp``` for predicting the impact of mutations on kinase-substrate phosphorylation
  - ```panoply_quilts``` for generating custom (personalized) protein databases
  - ```panoply_ptm_normalization``` for normalizing PTM (site-level) data to protein levels
* Updates to Multi-Omic NMF clustering (```panoply_mo_nmf```):
  - Balancing of -omes to prevent over-representation of one or more data types
  - Option for row/column z-scoring
  - Label-free data analysis
* RNA-protein correlation module ```panoply_rna_protein_correlation``` automatically row normalizes RNA data
* Consistent use of ```gene.id.col``` parameter to specify input data column with Gene Symbols
* Moved panda notebook to ```gcr.io``` to enable faster docker download
* panda notebook can run in Terra or via the command line
* Bug fixes, usability and robustness improvements


**v1.0**

Initial release
