#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
task panoply_ssgsea {

	File input_ds
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

    ## VM parameters
	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail
		
		# prepare GCT file
		/home/pgdac/src/preprocessGCT.R -i ${input_ds} -y ${yaml_file} -l ${default=NA level} -t ${default=NA id_type} -o ${default=NA id_type_out} -a ${default=NA acc_type} -s ${default=NA seqwin_col} --gene_symbol_column ${default=NA gene_col} -k ${default=NA humanize_gene}  -v ${default=NA SGT_col} -d ${default=NA loc} -m ${default=NA mode} -r "${default=NA mod_res}" -p '${default=NA mod_type}' -u TRUE -z /home/pgdac/src
	}

	output {
		# Outputs defined here
		File results=select_first(read_lines("fn.out")) # read processed GCT filename from fn.out
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
