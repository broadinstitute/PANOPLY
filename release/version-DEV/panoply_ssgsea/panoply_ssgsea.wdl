#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
task panoply_ssgsea {

	File input_ds
	File gene_set_database
	File yaml_file
	
        ## parameters to create gene-centric or single-site-centric 
        ## GCT files for ssGSEA / PTM-SEA
	String? level
 	String? id_type
	String? id_type_out
	String? acc_type
	String? seqwin_col
	String? gene_col
	Boolean? humanize_gene
	String? SGT_col
	Boolean? loc
	String? mode
	String? mod_res
	String? mod_type
	Boolean? preprocess_gct
	
	## ssGSEA / PTM-SEA parameters below	
	String output_prefix='results-ssgsea'
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
		
		# prepare GCT file
		/home/pgdac/src/preprocessGCT.R -i ${input_ds} -y ${yaml_file} -l ${default=NA level} -t ${default=NA id_type} -o ${default=NA id_type_out} -a ${default=NA acc_type} -s ${default=NA seqwin_col} --gene_symbol_column ${default=NA gene_col} -k ${default=NA humanize_gene}  -v ${default=NA SGT_col} -d ${default=NA loc} -m ${default=NA mode} -r "${default=NA mod_res}" -p '${default=NA mod_type}' -u ${default=NA preprocess_gct} -z /home/pgdac/src
		
		# update path to input_ds
		input_ds_proc=`cat fn.out`
		
		# run ssgsea/ptm-sea
		# don't use curly brackets for $input_ds_proc because it this is not a WDL variable 
		/home/pgdac/ssgsea-cli.R -i $input_ds_proc -y ${yaml_file} -d ${gene_set_database} -o ${default=NA output_prefix} -n ${default=NA sample_norm_type} -w ${default=NA weight} -c ${default=NA correl_type} -t ${default=NA statistic} -s ${default=NA output_score_type} -p ${default=NA nperm} -m ${default=NA min_overlap} -g ${default=NA global_fdr}

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
		docker : "broadcptacdev/panoply_ssgsea:DEV"
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
