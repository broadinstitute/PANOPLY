#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#

task panoply_metaboanalyst {
	File input_gct

	String output_prefix="results_metaboanalyst"
	File? yaml_file

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail
		
		Rscript /prot/proteomics/Projects/PGDAC/src/nmf.r -d ${sep="," ome_gcts} -o ${sep="," ome_labels} ${"-f " + sd_filt_min} ${"-g " + sd_filt_mode} ${"-u " + z_score} ${"-v " + z_score_mode} ${"-a " + gene_column} ${"-i " + organism_id} ${"--kmin " + kmin} ${"--kmax " + kmax} ${"-e " + exclude_2} ${"-n " + nrun} ${"-m " + nmf_method} ${"-s " + seed} -x ${output_prefix} ${"-y " + yaml_file} --libdir /prot/proteomics/Projects/PGDAC/src/
	}

	output {
		File results="${output_prefix}_NMF_results.tar.gz" # tar w/ expr GCT files + res.rank & parameters .Rdata files
		Int nclust=read_int("nmf_best_rank.txt")
		File? preprocess_figs="NMF_preprocessing_figures.tar.gz"
	}

	runtime {
		docker : "broadcptacdev/panoply_metaboanalyst:latest"
		memory : select_first ([memory, 32]) + "GB"
		disks : "local-disk " + select_first ([disk_space, 20]) + " HDD"
		cpu : select_first ([num_threads, 32]) + ""
		preemptible : select_first ([num_preemtions, 0])
	}

	meta {
		author : "C.M. Williams"
		email : "wcorinne@broadinstitute.org"
	}

}

################################################
## workflow
workflow panoply_metaboanalyst_workflow {
    call panoply_metaboanalyst
}
