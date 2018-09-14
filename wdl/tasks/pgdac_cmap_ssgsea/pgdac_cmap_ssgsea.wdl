task pgdac_cmap_ssgsea {
  # task adapted from pgdac_ssgsea; many inputs are set to specfic values for CMAP analysis
	File input_ds
	File gene_set_database
	Int? permutation_num
	String output_prefix = "${basename (input_ds, '.gctx')}" + "${if defined (permutation_num) then '-'+permutation_num else ''}"
	
	# other ssgsea options (below) are fixed for CMAP analysis
	String sample_norm_type = "rank"
	String correl_type = "rank"
	String statistic = "Kolmogorov-Smirnov"
	String output_score_type = "NES"

	Float weight = 0.0
	Int min_overlap = 5
	Int nperm = 0
	String global_fdr = "TRUE"
	String export_sigs = "FALSE"
	String ext_output = "FALSE"


  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions
	
	command {
		set -euo pipefail
		/home/pgdac/ssgsea-cli.R -i ${input_ds} -d ${gene_set_database} -o ${output_prefix} -n ${sample_norm_type} -w ${weight} -c ${correl_type} -t ${statistic} -s ${output_score_type} -p ${nperm} -m ${min_overlap} -g ${global_fdr} -e ${export_sigs} -x ${ext_output}
	}

	output {
		File scores="${output_prefix}-scores.gct"
	}

	runtime {
		docker : "broadcptac/pgdac_ssgsea:3"
    memory : select_first ([memory, 60]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, 16]) + ""
    preemptible : select_first ([num_preemptions, 2])
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}
}

workflow pgdac_cmap_ssgsea_workflow {
	call pgdac_cmap_ssgsea
}
