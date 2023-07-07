#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_internal_workflow/versions/1/plain-WDL/descriptor" as nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_sankey_workflow/versions/1/plain-WDL/descriptor" as sankey_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_assemble_results/versions/1/plain-WDL/descriptor" as assemble_wdl


################################################
##  workflow: nmf_balance (optional) + nmf + nmf_report + ssgsea + ssgsea_report
workflow panoply_nmf_workflow {
	String label

	# Data Upload
	Array[File]+ ome_gcts
	Array[String]+ ome_labels
	File? omes_tar # tar file with GCTs for analysis. trim filenames by default?

	# Auxilliary Files
	File gene_set_database # used for ssgsea
	File? yaml_file # used for colors
	File? groups_file # used for enrichement analysis
	
	# Analysis Parameters
	Float? tol
	Float? var
	String? zscore_mode

	# Module Toggles
	Boolean run_so_nmf
	Boolean run_mo_nmf
	Boolean run_sankey


	# zip ome_labels w/ ome files together
	Array[Pair[String,File]]+ ome_pairs = zip(ome_labels, ome_gcts)


	# Single-Ome NMF
	if (run_so_nmf) { # if we are running so_nmf
		scatter (pair in ome_pairs) {
			call nmf_wdl.panoply_nmf_internal_workflow as so_nmf { # call nmf_internal for each pair
				input:
					label="${label}-so_nmf",
					ome_pairs=pair,
					gene_set_database=gene_set_database,
					yaml_file=gene_set_database,
					groups_file=groups_file,

					balance_omes=false # no balancing for so_nmf
			}
		}
	}

	# Multi-Ome NMF
	if (run_mo_nmf) { # if we are running mo_nmf
		call nmf_wdl.panoply_nmf_internal_workflow as mo_nmf { # call nmf_internal for full array of pairs
			input:
				label="${label}-mo_nmf",
				ome_pairs=ome_pairs,
				gene_set_database=gene_set_database,
				yaml_file=gene_set_database,
				groups_file=groups_file
		}
	}

	# Sankey Diagrams
	if ( run_so_nmf && run_sankey ){ # if we ran so_nmf && want to run sankey
		call sankey_wdl.panoply_nmf_sankey_workflow as nmf_sankey {
			input:
				label = label,

				annot_files = so_nmf.membership,		# array of so-NMF results
				annot_file_labels = ,					# array of ome labels

				annot_file_primary = mo_nmf.membership, # single file with mo-NMF results
				annot_label_primary = "Multiomic",		# label for mo-NMF data

				annot_of_comparison="NMF.consensus",	# column for analysis
				annot_prefix="C"						# prefix to append to cluster values (e.g. 'C1' 'C2' 'C3')
		}
	}

	# Assemble Outputs into User-Friendly Structure
	call assemble_wdl.panoply_nmf_assemble_results as nmf_assemble {
		input:
			mo_nmf_tar = mo_nmf.nmf_tar,
			mo_nmf_report = mo_nmf.nmf_report,
			mo_nmf_ssgsea_tar = mo_nmf.nmf_ssgsea_tar,
			mo_nmf_ssgsea_report = mo_nmf.nmf_ssgsea_report,

			so_nmf_tar = so_nmf.nmf_tar,
			so_nmf_report = so_nmf.nmf_report,
			so_nmf_ssgsea_tar = so_nmf.nmf_ssgsea_tar,
			so_nmf_ssgsea_report = so_nmf.nmf_ssgsea_report

			sankey_tar = nmf_sankey.sankey_tar,
			sankey_report = nmf_sankey.sankey_report
	}


	output {
		File nmf_results=nmf_assemble.results
		File nmf_reports=nmf_assemble.report
		File? sankey_report=nmf_sankey.report
	}
}
