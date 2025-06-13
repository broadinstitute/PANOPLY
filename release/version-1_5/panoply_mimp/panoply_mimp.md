# ```panoply_mimp```

## Description

This module utilizes the MIMP algorithm for predicting the impact of missense mutations on kinase-substrate phosphorylation ([Wagih et al. 2015](https://www.nature.com/articles/nmeth.3396)). This module takes in mutation, phosphoproteomic, and protein sequence data and prepares them as input to the [rmimp R package](https://github.com/omarwagih/rmimp). For each sample individually, missense mutations in close proximity to phosphosites (called phosphorylation-related SNV, or "pSNV") are first identified, and then kinase specificity models are used to predict if each pSNV disrupts an existing phosphosite or creates a new phosphosite. These changes to expected phosphosites are referred to as *kinase rewiring events* as they can potentially alter known kinase networks. This module returns the results from individual samples as well as a summary of results across the whole input dataset, including heatmaps to visualize results. If a given sample does not have any MIMP results, this means that either there were no pSNVs identified or none of the pSNVs gave rise to kinase rewiring events.

## Input

### Required inputs:

* ```mutation_file```: (`.maf` file) mutation data table, must contain columns describing the mutation type (e.g. missense), the amino acid change and position, the sample identifier, and transcript identifier. 
* ```phospho_file```: (`.gct` file) phosphoproteomic data table, row metadata must include a columns for phosphosite and protein identifiers.
* ```fasta_file```: (`.fasta` file) list of protein sequences.
* ```ids_file```: (`.Rdata` file) mapping file between transcript and protein identifiers.
* ```master_yaml```: (`.yaml` file) master parameters file.
* ```output_prefix```: (String) prefix for naming the output tar file.

### Optional inputs:

* ```groups_file_path```: (`.csv` file, default = NULL) path to subset of sample annotations (groups file); if provided these annotations will be added as annotation tracks to output heatmaps.
* ```search_engine```: (String, default = "SpectrumMill") name of search engine used to process phosphoproteomic data. Options are "SpectrumMill" or "other." If "SpectrumMill," the user does not have to specify values for ```phosphosite_col``` and ```protein_id_col	``` below; if "other," the user must specify values for ```phosphosite_col``` and ```protein_id_col	```.
* ```phosphosite_col```: (String, default = NULL) if search_engine = "other," provide name of column in phospho row metadata that indicates phosphosite position.
* ```protein_id_col```: (String, default = NULL) if search_engine = "other," provide name of column in phospho row metadata that indicates protein identifier/accession number.
* ```mutation_AA_change_colname```: (String, default = "Protein_Change")  name of column in mutation MAF file that indicates amino acid change due to mutation.
* ```mutation_type_col```: (String, default =  "Variant_Classification")  name of column in mutation MAF file that indicates type of mutation (e.g. missense).
* ```sample_id_col```: (String, default =  "Tumor_Sample_Barcode") name of column in mutation MAF file that indicates sample identifier.
* ```transcript_id_col```: (String, default =  "refseq_mrna_id") name of column in mutation MAF file that indicates transcript identifier.

## Output

An output tar file called output_prefix_mimp_output.tar contains all results and summaries in a directory called "mimp_results_dir":

Results for individual samples are included in the "results_by_sample" sub-directory. Each sample has a further sub-directory labeled by sample identifier. These results include:
* In the "mutation_info" directory:
	* If the sample had mutations, full mutation data for the given sample (before any formatting): "sampleID_mutation_data_all.csv"
	* If any protein identifiers appeared in the phospho data that were not also in the fasta file, they were removed from analysis and listed here: "sampleID_phospho_identifiers_not_in_fasta.csv"
	* If there were any discrepancies between the reference amino acid in the mutation data and the amino acid at that position in the fasta file, they are listed here: "sampleID_mismatchedAA_mutation_vs_fasta.csv"
	* If there were any mutations that caused a gain of an S, T, or Y amino acid, they are listed here: "sampleID_STY_gain_mutation.csv"
	* If there were any mutations that caused a loss of an S, T, or Y amino acid, they are listed here: "sampleID_STY_loss_mutation.csv"
* In the "mimp_input" directory:
	* Formatted mutation data file used as input to the MIMP algorithm: "sampleID_mutation_mimp_input.csv"
	* Formatted phospho data file used as input to the MIMP algorithm: 01BR026_phospho_mimp_input.csv
* In the "mimp_results" directory:
	* If MIMP predicted any kinase rewiring events in the sample, full MIMP results for the sample: "sampleID_mimp_output_kinase_rewiring_events.csv"
	* If mimp predicted any kinase rewiring events in the sample, HTML report summarizing full MIMP results for the sample: "html/MIMP_results.html"


MIMP results for the whole dataset, if any, are included in the "mimp_results_dir" home directory:
* Concatenated list of full MIMP results from all samples: "mimp_output_kinase_rewiring_events_all.csv"
* Matrix of log-ratio results for all predicted kinase rewiring events by MIMP, indicating for each sample the likelihood and direction of altered kinase activity due to a specific point mutation: "kinase_rewiring_events_matrix_kinase_gene_mut_level.csv" (table); "mimp_results_predicted_kinase_rewiring_kinase_gene_mut_level_heatmap.pdf" and ".png" (heatmap)
* Matrix of log-ratio results summarizing overall altered kinase activity predicted by MIMP, indicating the maximum absolute value likelihood and direction of a kinase's altered activity in each sample: "kinase_rewiring_events_matrix_kinase_level.csv" (table); "mimp_results_predicted_kinase_rewiring_kinase_level_heatmap.pdf" and ".png" (heatmap)

Other files with mutation info for the full dataset in the "mimp_results_dir" home directory, if applicable:
* Concatenated list of all mutations causing a gain of an S, T, or Y amino acid in all samples: "all_STY_gain_mutation.csv"
* Concatenated list of all mutations causing a loss of an S, T, or Y amino acid in all samples: "all_STY_loss_mutation.csv"
* Concatenated list of all discrepancies between the reference amino acid in the mutation data and the amino acid at that position in the fasta file in all samples: "all_mismatchedAA_mutation_vs_fasta.csv"