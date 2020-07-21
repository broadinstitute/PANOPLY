<<<<<<< HEAD:hydrant/tasks/pgdac_mo_nmf/pgdac_mo_nmf.wdl
#########################################################################
## NMF
task pgdac_mo_nmf {

	#Inputs defined here
	File tar_file
	File yaml_file
	
	Int? kmin
	Int? kmax
	Int? nrun
	
	String? seed
	#String? cat_variable
	#String? cont_variable
	#String? class_colors

        String output_prefix="results_nmf"
	
	String? gene_column

	Boolean? no_plot
	Boolean? z_score
	Boolean? impute
	#Boolean? bnmf
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
		docker : "broadcptacdev/pgdac_mo_nmf:15"
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

#################################################################################
## ssGSEA
task pgdac_ssgsea {

	File input_ds
	
	## parameters to create gene-centric or single-site-centric 
	## GCT files for ssGSEA / PTM-SEA
	String level
	String? id_type
	String? id_type_out
	String? acc_type
	String? seqwin_col
	String? gene_col
	String? SGT_col
	Boolean? loc
	String? mode
	String? mod_res
	String? mod_type
	Boolean? preprocess_gct
	
	## ssGSEA / PTM-SEA parameters below
	File gene_set_database
	String output_prefix
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
		
		## output prefix
		##out_pref=${default="out" output_prefix}


		# prepare GCT file
		/home/pgdac/src/preprocessGCT.R -i ${input_ds} -l ${default="ssc" level} -t ${default="sm" id_type} -o ${default="seqwin" id_type_out} -a ${default="refseq" acc_type} --seqwin_column ${default="VMsiteFlanks" seqwin_col} --gene_symbol ${default="geneSymbol" gene_col} -v ${default="subgroupNum" SGT_col} -d ${default=true loc} -m ${default="median" mode} -r "${default="S|T|Y" mod_res}" -p ${default="p" mod_type} -u ${default=true preprocess_gct} -z /home/pgdac/src
	
		# update path to input_ds
		input_ds2=`cat fn.out`
		
		# default value for 'min.overlap'
		#min_overlap2=10
		#if [ ${default="ssc" level} = "ssc" ]; then
		#min_overlap2=5
		#fi

		# run ssgsea/ptm-sea       
		/home/pgdac/ssgsea-cli.R -i $input_ds2 -d ${gene_set_database} -o ${output_prefix} -n ${default="rank" sample_norm_type} -w ${default='0.75' weight} -c ${default="z.score" correl_type} -t ${default="area.under.RES" statistic} -s ${default="NES" output_score_type} -p ${default=1000 nperm} -m ${default=5 min_overlap} --globalfdr ${default=true global_fdr}
		
		# archive result files
		find * -regextype posix-extended -regex "^signature_gct/.*.gct$|^${output_prefix}.*.gct$|^.*.log.txt$|^.*parameters.txt$" -print0 | tar -czvf ${output_prefix}.tar.gz --null -T -
		}

	output {
		# Outputs defined here
		File results="${output_prefix}.tar.gz"
		}

	runtime {
		docker : "broadcptac/pgdac_ssgsea:5"
		memory : select_first([memory, 4]) + " GB"
		disks : "local-disk " + select_first([disk_space, 15]) + " HDD"
		cpu : select_first([num_threads, 2])
		preemptible : select_first ([num_preemtions, 2])
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}

}

################################################
##  workflow: mo-nmf plus ssgsea
workflow pgdac_mo_nmf_workflow {
	
	String? gene_col
	String level="gc"
	File gene_set_database
        File tar_file

	call pgdac_mo_nmf {
		input:
			tar_file=tar_file,
			gene_column=select_first([gene_col, "geneSymbol"])
	}

	call pgdac_ssgsea {
		input: 
			input_ds=pgdac_mo_nmf.feature_matrix_w,
			gene_set_database=gene_set_database,
			level=level,
			gene_col=select_first([gene_col, "geneSymbol"]),
			weight=1,
			mode="abs.max"
	}

	output {
		File nmf_clust=pgdac_mo_nmf.results
		File nmf_ssgsea=pgdac_ssgsea.results
	}
}
=======
task panoply_mo_nmf {
	#Inputs defined here
	File tar_file
	
	Int kmin
	Int kmax
	Int nrun
	
	String seed
	String class_variable
	String other_variable
	String class_colors
	String output_prefix
	String genes_of_interest
	String gene_column='geneSymbol'
  String lib_dir='/home/pgdac/src/'

	Boolean no_plot
	Boolean z_score
	Boolean impute
	
	Float sd_min=0

	Int memory
	Int disk_space
	Int num_threads
	Int num_preemtions
	
	command {
		set -euo pipefail
		#Command goes here
		Rscript /home/pgdac/src/mo-nmf.r -t ${tar_file} -l ${kmin} -m ${kmax} -n ${nrun} -s ${seed} -c ${class_variable} -o ${other_variable} -d ${class_colors} -r ${no_plot} -g ${genes_of_interest} -f ${sd_min} -b ${z_score} -i ${impute} -a ${gene_column} -z ${lib_dir}
		find * -type d -name "[0-9]*" -print0 | tar -czvf ${output_prefix}.tar --null -T -
		}

	output {
		#Outputs defined here
		File results="${output_prefix}.tar"
		}

	runtime {
		docker : "broadcptac/panoply_mo_nmf:5"
		memory : select_first ([memory, 7]) + "GB"
		disks : "local-disk " + select_first ([disk_space, 10]) + " SSD"
		cpu : select_first ([num_threads, 12]) + ""
		preemptible : select_first ([num_preemtions, 0])
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}

}

workflow panoply_mo_nmf_workflow {
	call panoply_mo_nmf
}
>>>>>>> dev:hydrant/tasks/panoply_mo_nmf/panoply_mo_nmf.wdl
