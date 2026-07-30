#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
version 1.1

task panoply_ssgsea {
	input {

		File input_ds
		File gene_set_database
		File yaml_file
		String output_prefix='results-ssgsea'

		## ssGSEA / PTM-SEA parameters below	
		String? sample_norm_type
		String? correl_type
		String? statistic
		String? output_score_type
		Float? weight
		Int? min_overlap
		String? tolerate_min_overlap_err # boolean value: should the WDL tolerate "not-enough-overlap" errors?
		Int? nperm
		Boolean? global_fdr

	    ## VM parameters
		Int? memory
		Int? disk_space
		Int? num_threads
		Int? num_preemptions
	}

	command {
		set -euo pipefail
		
		# run ssgsea/ptm-sea
		/home/pgdac/ssgsea-cli.R -i ${input_ds} -y ${yaml_file} -d ${gene_set_database} -o ${select_first([output_prefix, "NA"])} -n ${select_first([sample_norm_type, "NA"])} -w ${select_first([weight, "NA"])} -c ${select_first([correl_type, "NA"])} -t ${select_first([statistic, "NA"])} -s ${select_first([output_score_type, "NA"])} -p ${select_first([nperm, "NA"])} -m ${select_first([min_overlap, "NA"])} ${"-q " + tolerate_min_overlap_err} -g ${select_first([global_fdr, "NA"])} -z /home/pgdac


		## tar results

		# copy ${input_ds} to PWD, so it gets tarred with outputs
		cp ${input_ds} .

		# create regex to locate relevant outputs to tar
		result_regexpr="^signature_gct/.*.gct$|^${output_prefix}.*.gct$|^${basename(input_ds)}$|^.*.log.txt$|^.*parameters.txt$"
		find * -regextype posix-extended -regex $result_regexpr -print0 | tar -czvf ${output_prefix}.tar.gz --null -T -
		}

	output {
		# Outputs defined here
		File results="${output_prefix}.tar.gz"
		Boolean ssgsea_min_overlap_err=read_boolean("geneset_overlap_below_min.txt")
		}

	runtime {
		docker : "broadcptacdev/panoply_ssgsea:latest"
		memory : select_first ([memory, 8]) + "GB"
		disks : "local-disk " + select_first ([disk_space, 10]) + " HDD"
		cpu : select_first ([num_threads, 8]) + ""
		preemptible : select_first ([num_preemptions, 2])
		}

	meta {
		author : "Karsten Krug"
		email : "proteogenomics@broadinstitute.org"
	}

}

workflow panoply_ssgsea_workflow {
	call panoply_ssgsea

}