task pgdac_ssgsea {
	#Inputs defined here
	File input_ds
	File gene_set_database
	
	String output_prefix
	String sample_norm_type
	String correl_type
	String statistic
	String output_score_type

	Float weight
	Int min_overlap
	Int nperm
	Boolean global_fdr

	Int memory
	Int disk_space
	Int num_threads
	
	command {
		set -euo pipefail
		#Command goes here
		/home/pgdac/ssgsea-cli.R -i ${input_ds} -d ${gene_set_database} -o ${output_prefix} -n ${sample_norm_type} -w ${weight} -c ${correl_type} -t ${statistic} -s ${output_score_type} -p ${nperm} -m ${min_overlap} -g ${global_fdr}
		find * -regextype posix-extended -regex "^signature_gct/.*.gct$|^${output_prefix}.*.gct$|^.*.log.txt$|^.*parameters.txt$" -print0 | tar -czvf ${output_prefix}.tar.gz --null -T -
		}

	output {
		#Outputs defined here
		File results="${output_prefix}.tar.gz"
		}

	runtime {
		docker : "broadcptac/pgdac_ssgsea:4"
		memory : "${memory}GB"
		disks : "local-disk ${disk_space} HDD"
		cpu : "${num_threads}"
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}

}

workflow pgdac_ssgsea_workflow {
	call pgdac_ssgsea
}
