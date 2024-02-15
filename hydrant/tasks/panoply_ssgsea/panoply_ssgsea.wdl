#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
task panoply_ssgsea {

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
	Int? nperm
	Boolean? global_fdr

    ## VM parameters
	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail
		
		# run ssgsea/ptm-sea
		/home/pgdac/ssgsea-cli.R -i ${input_ds} -y ${yaml_file} -d ${gene_set_database} -o ${default=NA output_prefix} -n ${default=NA sample_norm_type} -w ${default=NA weight} -c ${default=NA correl_type} -t ${default=NA statistic} -s ${default=NA output_score_type} -p ${default=NA nperm} -m ${default=NA min_overlap} -g ${default=NA global_fdr} -z /home/pgdac

		# set wdl variable 'output_prefix' to the value specified in the yaml file,
		# if not specified via cmd line 
		#output_prefix=${default="NA" output_prefix}
		output_prefix_local="results"
		#if( $output_prefix_local = "NA" ); then
 		#    output_prefix_local=`cat $yaml_file | grep output_prefix | awk '\\123print $2\\125'`
		#fi

		# archive result files
		result_regexpr="^signature_gct/.*.gct$|^${output_prefix}.*.gct$|^.*.log.txt$|^.*parameters.txt$"
		result_regexpr="$result_regexpr|^$input_ds_proc"
		
		find * -regextype posix-extended -regex $result_regexpr -print0 | tar -czvf ${output_prefix}.tar.gz --null -T -
		}

	output {
		# Outputs defined here
		File results="${output_prefix}.tar.gz"
		}

	runtime {
		docker : "broadcptacdev/panoply_ssgsea:latest"
		memory : select_first ([memory, 8]) + "GB"
		disks : "local-disk " + select_first ([disk_space, 10]) + " HDD"
		cpu : select_first ([num_threads, 8]) + ""
		preemptible : select_first ([num_preemtions, 2])
		}

	meta {
		author : "Karsten Krug"
		email : "proteogenomics@broadinstitute.org"
	}

}

workflow panoply_ssgsea_workflow {
	call panoply_ssgsea
}
