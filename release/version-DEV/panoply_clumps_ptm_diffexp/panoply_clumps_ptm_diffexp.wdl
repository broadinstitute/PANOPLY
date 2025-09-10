#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#


task panoply_clumps_ptm_diffexp {
	# PTM GCT files; must include at least one
	File? pSTY_gct
	File? acK_gct
	File? ubK_gct

	File groupsFile
	String? sample_id_col
  	Float? fdr_cutoff
	Int? min_samples
	Int? max_annot_levels

	String output_prefix="results"
	File yaml_file

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemptions
	
	command {
		set -euo pipefail

		Rscript /prot/proteomics/Projects/PGDAC/src/clumps_diffexp.r ${"-p " + pSTY_gct} ${"-a " + acK_gct} ${"-u " + ubK_gct} \
		-g ${groupsFile} ${"-s " + sample_id_col} ${"-f " + fdr_cutoff} ${"-m " + min_samples} ${"-l " + max_annot_levels} \
		-x ${output_prefix} ${"-y " + yaml_file} --libdir /prot/proteomics/Projects/PGDAC/src/ 
	}

	output {
		Array[File] diff_exp_files = glob( "${output_prefix}_*_diff_exp.tsv" ) # all differential-expression files
		File log_file = "${output_prefix}_log_file.csv"					 # log file listing which annotations were analyzed
	}

	runtime {
		docker : "broadcptacdev/panoply_clumps_ptm_diffexp:DEV"
		memory : select_first ([memory, 32]) + "GB"
		disks : "local-disk  " + select_first ([disk_space, 20]) + " HDD"
		cpu : select_first ([num_threads, 32]) + ""
		preemptible : select_first ([num_preemptions, 0])
	}

	meta {
		author : "C.M. Williams"
		email : "proteogenomics@broadinstitute.org"
	}

}

################################################
## workflow
workflow panoply_clumps_ptm_diffexp_workflow {
	call panoply_clumps_ptm_diffexp {
	}
}
