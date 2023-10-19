#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#

task panoply_nmf_postprocess {
    File nmf_results # rdata object with results of NMF()
    Int nclust # best number of cluster
    File expr_comb # expression matrix GCT file
    File expr_comb_nn # non-negative expression matrix GCT file

	String output_prefix="results_nmf"
	File? yaml_file # mostly for colors
	File? groups_file # for enrichment analysis

	String? gene_column
	Float? pval_signif

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail
		
		Rscript /prot/proteomics/Projects/PGDAC/src/nmf_postprocess.R --nmf_results ${nmf_results} --rank_top ${nclust} --expr_comb ${expr_comb} --expr_comb_nn ${expr_comb_nn} ${"-g " + groups_file} ${"-a " + gene_column} ${"-p " + pval_signif} -x ${output_prefix} ${"-y " + yaml_file} --libdir /prot/proteomics/Projects/PGDAC/src/

		tar -czvf nmf_results.tar ${output_prefix}_K${nclust}_* # tar all files with 
	}

	output {
		File results="nmf_results.tar"
		File membership="${output_prefix}_K${nclust}_clusterMembership.csv"
	}

	runtime {
		docker : "broadcptacdev/panoply_nmf_postprocess:latest"
		memory : select_first ([memory, 64]) + "GB"
		disks : "local-disk " + select_first ([disk_space, 20]) + " HDD"
		cpu : select_first ([num_threads, 32]) + ""
		preemptible : select_first ([num_preemtions, 0])
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}

}

################################################
workflow panoply_nmf_postprocess_workflow {
	call panoply_nmf_postprocess
}
