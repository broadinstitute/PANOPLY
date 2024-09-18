#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#

task panoply_nmf_postprocess {
    File nmf_results		# tar w/ expr GCT files + res.rank & parameters .Rdata files
    Int nclust				# best number of cluster

	String? gene_column
	File? groups_file		# for enrichment analysis

	Float? pval_signif		# significant p-value for cluster-enrichement
	Float? feature_fdr		# fdr threshold for driver-feature t-test
	Int? max_annot_levels	# max number of annotation-levels to allow for discrete variables
	Int? top_n_features		# max number of driver features (per cluster) to create expression-boxplots for 

	String output_prefix="output"
	File? yaml_file

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail
		
		Rscript /prot/proteomics/Projects/PGDAC/src/nmf_postprocess.R --nmf_results ${nmf_results} --rank_top ${nclust} ${"-g " + groups_file} ${"-a " + gene_column} ${"-p " + pval_signif} ${"-q " + feature_fdr} ${"-l " + max_annot_levels} ${"-t " + top_n_features} -x ${output_prefix} ${"-y " + yaml_file} --libdir /prot/proteomics/Projects/PGDAC/src/

	}

	output {
		File results="${output_prefix}_NMF_postprocess.tar.gz" # contains all outputs of NMF post-processing + parameters
		File membership="${output_prefix}_K${nclust}_clusterMembership.tsv"
		File feature_matrix_w=select_first(glob("${output_prefix}_K${nclust}_W_rowNorm_combined_signed_n*.gct")) # select first and only match (of array-length-1)
		Boolean ssgsea_viable=read_boolean("gene_col_toggle.txt")
	}

	runtime {
		docker : "broadcptacdev/panoply_nmf_postprocess:latest"
		memory : select_first ([memory, 64]) + "GB"
		disks : "local-disk " + select_first ([disk_space, 20]) + " HDD"
		cpu : select_first ([num_threads, 32]) + ""
		preemptible : select_first ([num_preemtions, 0])
	}

	meta {
		author : "C.M. Williams"
		email : "proteogenomics@broadinstitute.org"
	}

}

################################################
workflow panoply_nmf_postprocess_workflow {
	call panoply_nmf_postprocess
}
