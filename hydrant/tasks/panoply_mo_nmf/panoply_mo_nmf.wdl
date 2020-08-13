import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea/versions/9/plain-WDL/descriptor" as panoply_ssgsea_wdl

#########################################################################
## NMF
task panoply_mo_nmf {

	#Inputs defined here
	File tar_file
	File yaml_file
	
	Int? kmin
	Int? kmax
	Int? nrun
	
	String? seed

        String output_prefix="results_nmf"
	
	String? gene_column

	Boolean? no_plot
	Boolean? z_score
	Boolean? impute

	Boolean? exclude_2
	Float? sd_min

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail
		#Command goes here
		Rscript /home/pgdac/src/mo-nmf.r -t ${tar_file} -y ${yaml_file} -l ${default=NA kmin} -m ${default=NA kmax} -n ${default=NA nrun} -s ${default=NA seed} -r ${default=false no_plot} -f ${default=NA sd_min} -u ${default=NA z_score} -i ${default=NA impute} -a ${default=NA gene_column} -e ${default=NA exclude_2} -z /home/pgdac/src/
		
		find * -type d -name "[0-9]*" -print0 | tar -czvf ${output_prefix}.tar --null -T -
		
		find . -name "matrix_W_combined_signed.gct" | xargs cp -t .
		}

	output {
		#Outputs defined here
		File results="${output_prefix}.tar"
		File feature_matrix_w="matrix_W_combined_signed.gct"
		}

	runtime {
		docker : "broadcptacdev/panoply_mo_nmf:latest"
		memory : select_first ([memory, 8]) + "GB"
		disks : "local-disk " + select_first ([disk_space, 10]) + " HDD"
		cpu : select_first ([num_threads, 10]) + ""
		preemptible : select_first ([num_preemtions, 2])
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}

}

################################################
##  workflow: mo-nmf plus ssgsea
workflow panoply_mo_nmf_workflow {
	
	File tar_file
	File gene_set_database
        File yaml_file

	call panoply_mo_nmf {
		input:
			tar_file=tar_file,
			yaml_file=yaml_file
	}

	call panoply_ssgsea_wdl.panoply_ssgsea {
		input:
			input_ds=panoply_mo_nmf.feature_matrix_w,
			gene_set_database=gene_set_database,
			yaml_file=yaml_file,
			output_prefix='nmf-ssgsea',
			weight=1
	}

	output {
		File nmf_clust=panoply_mo_nmf.results
		File nmf_ssgsea=panoply_ssgsea.results
	}
}
