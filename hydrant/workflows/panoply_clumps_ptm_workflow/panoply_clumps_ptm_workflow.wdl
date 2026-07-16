#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
version 1.0

import "../../tasks/panoply_clumps_ptm_diffexp/panoply_clumps_ptm_diffexp.wdl" as diffexp_wdl
import "../../tasks/panoply_clumps_ptm_mapping/panoply_clumps_ptm_mapping.wdl" as mapping_wdl
import "../../tasks/panoply_clumps_ptm/panoply_clumps_ptm.wdl" as analysis_wdl
import "../../tasks/panoply_clumps_ptm_postprocess/panoply_clumps_ptm_postprocess.wdl" as postprocess_wdl
import "../../tasks/panoply_clumps_ptm_report/panoply_clumps_ptm_report.wdl" as report_wdl

################################################
##  workflow: panoply_clumps_ptm_diffexp + panoply_clumps_ptm_mapping + panoply_clumps_ptm
workflow panoply_clumps_ptm_workflow {
	# PTM GCT files; must include at least one
	input {
		File? pSTY_gct
		File? acK_gct
		File? ubK_gct

		File groupsFile

		## Default Database Files		
		# Google-Cloud Bucket with PDB Directory, split into tarfiles
		String PDB_ref_bucket = "gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/pdbs_2025_tars/"
		# Uniprot FASTA reference file
		File UNIPROT_SWISSPROT = "gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/reference_files/uniprot_sprot.fasta"
		# SIFTS database
		File SIFTS_DB = "gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/reference_files/pdb_chain_uniprot.tsv"

		File? FASTA_ref_file

		String? accession_col
		String? variable_sites_col

		File? mapping_file					# pre-generated mapping file, to skip mapping module
		File? mapping_params				# parameter file from mapping

		String output_prefix
		File yaml_file
	}

	call diffexp_wdl.panoply_clumps_ptm_diffexp as diffexp {
		input:
			pSTY_gct = pSTY_gct,
			acK_gct = acK_gct,
			ubK_gct = ubK_gct,
			groupsFile = groupsFile,
			yaml_file = yaml_file,
			output_prefix = output_prefix
	}

    if ( !defined(mapping_file) ){ # only run if a mapping file is not provided
		call mapping_wdl.panoply_clumps_ptm_mapping as mapping {
			input:
				pSTY_gct = pSTY_gct,
				acK_gct = acK_gct,
				ubK_gct = ubK_gct,
				PDB_ref_bucket = PDB_ref_bucket,
				UNIPROT_SWISSPROT = UNIPROT_SWISSPROT,
				SIFTS_DB = SIFTS_DB,
				FASTA_ref_file = select_first([FASTA_ref_file]),
				accession_col = accession_col,
				variable_sites_col = variable_sites_col,
				yaml_file = yaml_file,
				output_prefix = output_prefix
		}
	}

  	scatter (diff_exp in diffexp.diff_exp_files) {
		call analysis_wdl.panoply_clumps_ptm as analysis {
		    input:
		    	diff_exp_file = diff_exp,
		    	var_sites_file = select_first([mapping_file, mapping.filt_results]),
		        PDB_ref_bucket = PDB_ref_bucket,
				accession_col = accession_col,
				variable_sites_col = variable_sites_col,
	        	yaml_file = yaml_file,
		        output_prefix = sub(basename(diff_exp), "_diff_exp\\.tsv$", "") # use diffexp file as output prefix (to keep annotation)
		}

		call postprocess_wdl.panoply_clumps_ptm_postprocess as postprocess {
		    input:
		    	results_tar = analysis.results,
	        	yaml_file = yaml_file,
		        output_prefix = sub(basename(diff_exp), "_diff_exp\\.tsv$", "") # use diffexp file as output prefix (to keep annotation)
		}
	}

	call report_wdl.panoply_clumps_ptm_report as report {
	    input:
	    	postprocess_results = postprocess.results, 			 # array of results AND figures from postprocess module
	    	mapping_params = select_first([mapping_params, mapping.mapping_params]),
	    	label = output_prefix
	}

	output{
	    Array[File] clumps_ptm_results = postprocess.results     # array of results AND figures from postprocess module
	    File clumps_ptm_report = report.report
	}
}
