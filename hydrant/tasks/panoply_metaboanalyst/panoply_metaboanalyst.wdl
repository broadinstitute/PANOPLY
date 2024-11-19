#
# Copyright (c) 2024 The Broad Institute, Inc. All rights reserved.
#

task panoply_metaboanalyst {
	File meta_gct
	String? meta_id_type
	File? omic_gct
	String? ome_type
	String? gene_column
	String? gene_id_type

	String? anal_type
	String? pval_comb

	Int? max_annot_levels
	String? pval_signif
	Int? top_n_networks

	String output_prefix="results_metaboanalyst"
	File? groups_file
	File yaml_file

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail

		Rscript /prot/proteomics/Projects/PGDAC/src/MetaboAnalyst_script.R '--metabolome_gct' ${meta_gct} ${'--meta_id_type ' + meta_id_type} ${'--ome_gct ' + omic_gct} ${'--ome_type ' + ome_type} ${'--gene_column ' + gene_column} ${'--gene_id_type ' + gene_id_type} ${'--groups_file ' + groups_file} ${"--max_annot_levels " + max_annot_levels} ${"--anal_type " + anal_type} ${"--pval_comb " + pval_comb} ${"--pval_signif " + pval_signif} ${"--top_n_networks " + top_n_networks} ${"--output_prefix " + output_prefix} ${"--yaml " + yaml_file} --libdir /prot/proteomics/Projects/PGDAC/src/
	}

	output {
		File results="${output_prefix}_MetaboAnalyst.tar.gz" # tar w/ outut files
	}

	runtime {
		docker : "broadcptacdev/panoply_metaboanalyst:latest"
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
workflow panoply_metaboanalyst_workflow {
    call panoply_metaboanalyst
}
