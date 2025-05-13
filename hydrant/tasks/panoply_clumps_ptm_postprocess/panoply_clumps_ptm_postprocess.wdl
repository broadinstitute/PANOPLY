#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#


task panoply_clumps_ptm_postprocess {
	File results_tar							# tar-file with results from panoply_clumps_ptm

	Float? fdr_threshold

	String output_prefix="results"
	File yaml_file

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail

		python -u /prot/proteomics/Projects/PGDAC/src/clumps_ptm_postprocess.py \
			--results_tar ${results_tar} -y ${yaml_file} --output_prefix ${output_prefix} \
			${'--fdr_threshold ' + fdr_threshold}

		# tar full results
		tar -czf ${output_prefix}_clumps_ptm_full_results.tar -C ${output_prefix}/ . # tar results
		# tar just the figures
		tar -czf ${output_prefix}_clumps_figures.tar -C ${output_prefix}/figures . # tar results
	}

	output {
		File results="${output_prefix}_clumps_ptm_full_results.tar"
		File figures_only="${output_prefix}_clumps_figures.tar"
	}

	runtime {
		docker : "broadcptacdev/panoply_clumps_ptm_postprocess:latest"
		memory : select_first ([memory, 32]) + "GB"
		disks : "local-disk  " + select_first ([disk_space, 20]) + " HDD"
		cpu : select_first ([num_threads, 32]) + ""
		preemptible : select_first ([num_preemtions, 0])
	}

	meta {
		author : "C.M. Williams"
		email : "proteogenomics@broadinstitute.org"
	}

}

################################################
## workflow
workflow panoply_clumps_ptm_postprocess_workflow {
	call panoply_clumps_ptm_postprocess {
	}
}
