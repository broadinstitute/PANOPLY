#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#

task panoply_nmf {
    Array[File]+ ome_gcts
    Array[String]+ ome_labels

	String output_prefix="results_nmf"
	File? yaml_file

	Float? sd_filt_min
	String? sd_filt_mode
	String? z_score_mode
	String? gene_column
	
	Int? kmin
	Int? kmax
	Int? nrun
	# Boolean? bayesian
	String? nmf_method # in the YAML
	String? seed

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail
		
		Rscript /prot/proteomics/Projects/PGDAC/src/nmf.r -d ${sep="," ome_gcts} -o ${sep="," ome_labels} ${"-f " + sd_filt_min} ${"-g " + sd_filt_mode} ${"-v " + z_score_mode} ${"-a " + gene_column}  ${"--kmin " + kmin} ${"--kmax " + kmax} ${"-n " + nrun} ${"-m " + nmf_method} ${"-s " + seed} -x ${output_prefix} ${"-y " + yaml_file} --libdir /prot/proteomics/Projects/PGDAC/src/
	}

	output {
		File results="nmf_res.Rdata"
		File gct_comb=select_first(glob("${output_prefix}_combined_n*.gct")) # select first/only match of array-length-1
		File gct_comb_nn=select_first(glob("${output_prefix}_combinedNonNegative*.gct")) # select first/only match of array-length-1
		Int nclust=read_int("nmf_best_rank.txt")
	}

	runtime {
		docker : "broadcptacdev/panoply_nmf:latest"
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
##  workflow: mo-nmf
workflow panoply_nmf_workflow {
	call panoply_nmf
}
