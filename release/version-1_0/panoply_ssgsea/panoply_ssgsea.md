Documentation at https://github.com/broadinstitute/PANOPLY/blob/release-1_0/release/version-1_0/panoply_ssgsea/panoply_ssgsea.md

# ```panoply_ssgsea```

## Description

Performs single sample Gene Set Enrichment Analysis (ssGSEA) or PTM-Signature Enrichment Analysis (PTM-SEA) [1] on each column of the input data matrix. This module is based on the implementation available at the [ssGSEA2.0 GitHub repository](https://github.com/broadinstitute/ssGSEA2.0). 

This is an updated version of the original ssGSEA [2,3] R-implementation. Depending on the input dataset and chosen database (gene sets or PTM signatures), the software performs either ssGSEA or PTM-SEA, respectively. The Molecular Signatures Database ([MSigDB](http://software.broadinstitute.org/gsea/msigdb/)) [4] provides a large collection of curated gene sets.  Gene sets are stored as plain text in  [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) format. A current version of MSigDB gene set collections can be found in the ```db/msigdb``` subfolder. MSigDB gene sets are realeased under [Creative Commons Attribution 4.0 International License](http://software.broadinstitute.org/gsea/msigdb_license_terms.jsp). The license terms can be found in the```db/msigdb``` folder.

File formats supported by ssGSEA2.0/PTM-SEA are Gene Cluster Text [GCT v1.2](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29) or [GCT v1.3](https://clue.io/connectopedia/gct_format) files. [Morpheus](https://software.broadinstitute.org/morpheus/) provides a convenient way to convert your data tables into GCT format.

For more information about the GSEA method and MSigDB please visit http://software.broadinstitute.org/gsea/.



## Input

### Required inputs:

* ```input_ds```: (`.gct` file) input GCT file
* ```gene_set_database```: (`.gmt` file) gene set database
* ```yaml_file```: (`.yaml` file) master-parameters.yaml

### Optional inputs:

#### preprocessing
* ```preprocess_gct```: (Boolean) If FALSE preprocessing will be skipped and the GCT file will be used as is (default: FALSE).
* ```acc_type```: (String) Type of accession number in 'rid' object in GCT file ("uniprot", "refseq" (default), "symbol").
* ```id_type```: (String) Notation of site-ids: 'sm' - Spectrum Mill (default); 'wg' - Web Gestalt; 'ph' - Philosopher. Only relevant for PTM-SEA.
* ```id_type_out``` (String) Type of site id for output: 'uniprot'(default), 'refseq', 'seqwin'. Only relevant for PTM-SEA.
* ```level``` (String) Mode of report:
  + 'ssc' - single-site-centric
  + 'gc' - gene-centric (default)
  + 'gcr' - gene-centric-redundant
* ```loc``` (Boolean) If TRUE only fully localized sites will be considered (default: TRUE). Localization infromation is expected to be encoded in the site identifier. Respective parsing rules are determined by '--id_type'.
* ```gene_col```: (String) Name of column listing gene names; used for gene centric reports (default: "geneSymbol").
* ```seqwin_col```: (String) "Column containing flanking sequences, separated by '|'. Only relevant for PTM-SEA and if '--id_type_out' = 'seqwin' (default: 'VMsiteFlanks').
* ```SGT_col```: (String) Column used to collpase subgroup-top (SGT) reports (default: "subgroupNum). Only relevant for Spectrum Mill protein reports.
* ```mod_res```: (String) Modified residues, e.g. "S|T|Y" or "K" (default: "S|T|Y").
* ```mod_type```: (String) Type of post-translational modification, e.g "p" for phospho (default) or "ac" for acetylation
* ```mode```: (String) Determines how multiple features (e.g. proteins, PTM sites, etc.) mapping to the same gene symbol will be aggregated: 
  + "mean" - mean
  + "median" - median
  + "sd - most variable (standard deviation) across sample columns
  + "SGT" - subgroup top: first subgroup in protein group (Spectrum Mill)
  + "abs.max" - for log-transformed, signed p-values"

#### ssGSEA/PTM-SEA
* ```correl_type```: (String) Correlation type: "z.score" (default), "rank", "symm.rank".
* ```global_fdr```: (Boolean) If TRUE global FDR across all data columns is calculated (default: FALSE).
* ```min_overlap```: (Integer) Minimal overlap between signature and data set (default: 10).
* ```nperm```: (Integer) Number of permutations (default: 1000).
* ```output_score_type```: (String) Score type: "ES" - enrichment score,  "NES" - normalized ES (default).
* ```statistic```: (String) Test statistic: "area.under.RES" (default), "Kolmogorov-Smirnov"
* ```weight```: (Float) When weight=0, all genes have the same weight; if weight>0 actual values matter and can change the resulting score (default: 0.75).
* ```sample_norm_type```: (String) Sample normalization: "rank"(default), "log", "log.rank"
* ```output_prefix```: (String) File prefix for output files.

## Output

* ```results```: (`.tar.gz` file) tarball including:

* `${output_prefix}-scores.gct`: GCT file with enrichment scores as data matrix (@mat)
* `${output_prefix}-pvalues.gct`: GCT file with nominal p-values as data matrix (@mat)
* `${output_prefix}-fdr-pvalues.gct`: GCT file with FDR-corrected p-values as data matrix (@mat)
* `${output_prefix}-combined.gct`: GCT file with enrichment scores as data matrix (@mat) as well as nominal and FDR-corrected p-values as row-description metadata (@rdesc).
* `${output_prefix}-parameters.txt`: Text file summarizing parameters.
* `${output_prefix}-ssgse.log.txt`: Text file tracking progess.
* `signature_gct`: Directory with individual GCT files (one for each gene set) with the original data used as input for ssGSEA/PTM-SEA as data matrix (@mat).

	
## References
1.  Krug, K., Mertins, P., Zhang, B., Hornbeck, P., Raju, R., Ahmad, R., . Szucs, M., Mundt, F., Forestier, D., Jane-Valbuena, J., Keshishian, H., Gillette, M. A., Tamayo, P., Mesirov, J. P., Jaffe, J. D., Carr, S. A., Mani, D. R. (2019). **A curated resource for phosphosite-specific signature analysis**, Molecular & Cellular Proteomics (in Press). http://doi.org/10.1074/mcp.TIR118.000943

1.  Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Susan, E., Dunn, I. F., . Hahn, W. C. (2010). **Systematic RNA interference reveals that oncogenic KRAS- driven cancers require TBK1**, Nature, 462(7269), 108-112. https://doi.org/10.1038/nature08460

1. Abazeed, M. E., Adams, D. J., Hurov, K. E., Tamayo, P., Creighton, C. J., Sonkin, D., et al. (2013).
       **Integrative Radiogenomic Profiling of Squamous Cell Lung Cancer.** Cancer Research, 73(20), 6289-6298.
       http://doi.org/10.1158/0008-5472.CAN-13-1616

1. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005).
   **Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.**
  Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. http://doi.org/10.1073/pnas.0506580102
	