#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
version 1.1

task panoply_preprocess_gct {
	input {

		File input_ds
		File yaml_file
		String output_prefix='input-ssgsea'

	    ## parameters to create gene-centric or single-site-centric
	    ## GCT files for ssGSEA / PTM-SEA
	    ## Defaults are read from the 'panoply_preprocess_gct' section of yaml_file;
	    ## any parameter set here overrides the yaml value.
		String? level
	 	String? id_type
		String? id_type_out
		String? acc_type
		String? organism # human, mouse, rat or zebrafish; required to map RefSeq / gene symbols to UniProt
		String? seqwin_col
		String? gene_col
		Boolean? humanize_gene
		String? SGT_col
		Boolean? loc
		String? mode
		String? mod_res
		String? mod_type

	    ## VM parameters
		Int? memory
		Int? disk_space
		Int? num_threads
		Int? num_preemptions
	}

	command {
		set -euo pipefail
		
		# prepare GCT file with ssGSEA3::preprocessGCT
		# (parameters that are not set are omitted, so preprocess-gct-cli.R takes them from yaml_file)
		Rscript /home/pgdac/src/preprocess-gct-cli.R \
			--input ${input_ds} --yaml ${yaml_file} \
			${"--level " + level} \
			${"--id_type " + id_type} \
			${"--id_type_out " + id_type_out} \
			${"--acc_type_in " + acc_type} \
			${"--organism " + organism} \
			${"--seqwin_column '" + seqwin_col + "'"} \
			${"--gene_symbol_column '" + gene_col + "'"} \
			${"--humanize_gene_names " + humanize_gene} \
			${"--sgt_column '" + SGT_col + "'"} \
			${"--localized " + loc} \
			${"--mode " + mode} \
			${"--residue '" + mod_res + "'"} \
			${"--ptm " + mod_type}

		mkdir file_output
		mv `cat fn.out` file_output/ # move the output file from preprocess-gct-cli.R into a folder
	}

	output {
		# Outputs defined here
		File result=glob('file_output/*.gct')[0]
	}

	runtime {
		docker : "broadcptacdev/panoply_preprocess_gct:test"
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

workflow panoply_preprocess_gct_workflow {
	call panoply_preprocess_gct

}
