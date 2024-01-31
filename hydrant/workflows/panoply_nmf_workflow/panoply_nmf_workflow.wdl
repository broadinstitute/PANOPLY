#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_select_all_pairs/versions/1/plain-WDL/descriptor" as select_pairs
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_internal_workflow/versions/12/plain-WDL/descriptor" as nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_sankey_workflow/versions/9/plain-WDL/descriptor" as sankey_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_assemble_results/versions/12/plain-WDL/descriptor" as assemble_wdl


################################################
##  workflow: mo_nmf + so_nmf + sankey + assemble
workflow panoply_nmf_workflow {
	String label

	## Data Upload
    Array[Pair[String?,File?]]+ ome_pairs_input
	# Array[File]+ ome_gcts		# array of GCT files
	# Array[String]+ ome_labels	# array of ome-labels for those GCT files (MUST MATCH ORDER)
	# File? omes_tar			# tar file with GCTs for analysis. not set up.

	## Auxilliary Files
	File gene_set_database		# used for ssgsea
	File yaml_file				# default parameters & figure colors
	File? groups_file			# datatable with annotations-of-interest (for figures & enrichement analysis)
	
	## Balance Parameters
	Float? tol
	Float? var
	String? zscore_mode
	
	## NMF Parameters
	Int? kmin
	Int? kmax
	String? exclude_2		# true / false
	String? nmf_method		# options in the YAML
	String? seed			# 'random' for random seed, or numeric for explicit seed

	## Module Toggles
	Boolean run_so_nmf
	Boolean run_mo_nmf
	Boolean run_sankey
    
	# select extant pairs from ome_pairs_input
    call select_pairs.panoply_select_all_pairs as select_pairs {
    	input:
        	pairs_input = ome_pairs_input
    }

	# Single-Ome NMF
	if (run_so_nmf) { # if we are running so_nmf
		scatter (pair in select_pairs.pairs) {
			call nmf_wdl.panoply_nmf_internal_workflow as so_nmf { # call nmf_internal for each pair
				input:
					label="${label}-so_nmf-${pair.left}",
					ome_labels=[pair.left],
					ome_gcts=[pair.right],
					gene_set_database=gene_set_database,
					yaml_file=yaml_file,
					groups_file=groups_file,
					kmin=kmin,
					kmax=kmax,
					exclude_2=exclude_2,
					nmf_method=nmf_method,
					seed=seed
			}
		}
	}

	# Multi-Ome NMF
	if (run_mo_nmf) { # if we are running mo_nmf
    	
		call nmf_wdl.panoply_nmf_internal_workflow as mo_nmf { # call nmf_internal for full array of pairs
			input:
				label="${label}-mo_nmf",
				ome_labels=select_pairs.pair_string,
				ome_gcts=select_pairs.pair_file,
				gene_set_database=gene_set_database,
				yaml_file=yaml_file,
				groups_file=groups_file,
				kmin=kmin,
				kmax=kmax,
				exclude_2=exclude_2,
				nmf_method=nmf_method,
				seed=seed
		}
	}

	# Sankey Diagrams
	if ( run_so_nmf && run_sankey ){ # if we ran so_nmf && want to run sankey
		call sankey_wdl.panoply_sankey_workflow as nmf_sankey {
			input:
				label = label,

				annot_files = so_nmf.nmf_membership,			## array of so-NMF results
				annot_file_labels = select_pairs.pair_string,	## array of ome labels

				annot_file_primary = mo_nmf.nmf_membership, 	## single file with mo-NMF results
				annot_label_primary = "Multiomic",				## label for mo-NMF data

				annot_column="NMF.consensus"					## column for analysis
		}
	}

	# Assemble Outputs into User-Friendly Structure
	call assemble_wdl.panoply_nmf_assemble_results as nmf_assemble {
		input:
			mo_nmf_results = mo_nmf.nmf_results,
			mo_nmf_figures = mo_nmf.nmf_figures,
			mo_nmf_report = mo_nmf.nmf_report,
			mo_nmf_ssgsea_tar = mo_nmf.nmf_ssgsea_tar,
			mo_nmf_ssgsea_report = mo_nmf.nmf_ssgsea_report,

			so_nmf_results = so_nmf.nmf_results,
			so_nmf_figures = so_nmf.nmf_figures,
			so_nmf_report = so_nmf.nmf_report,
			so_nmf_ssgsea_tar = so_nmf.nmf_ssgsea_tar,
			so_nmf_ssgsea_report = so_nmf.nmf_ssgsea_report,

			sankey_tar = nmf_sankey.sankey_tar,
			sankey_report = nmf_sankey.sankey_report_file,

			output_results_tar = "${label}_nmf_results.tar.gz",
			output_reports_tar = "${label}_nmf_reports.tar.gz"
	}


	output {
		File nmf_results=nmf_assemble.nmf_results
		File nmf_reports=nmf_assemble.nmf_reports
		# File? sankey_report=nmf_sankey.sankey_report_file
	}
}
