#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_internal_workflow/versions/10/plain-WDL/descriptor" as nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_sankey_workflow/versions/1/plain-WDL/descriptor" as sankey_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_assemble_results/versions/4/plain-WDL/descriptor" as assemble_wdl


################################################
##  workflow: mo_nmf + so_nmf + sankey + assemble
workflow panoply_nmf_workflow {
	String label

	## Data Upload
	Array[File]+ ome_gcts		# array of GCT files
	Array[String]+ ome_labels	# array of ome-labels for those GCT files (MUST MATCH ORDER)
	# File? omes_tar			# tar file with GCTs for analysis. not set up.

	## Auxilliary Files
	File gene_set_database		# used for ssgsea
	File yaml_file				# default parameters & figure colors
	File? groups_file			# datatable with annotations-of-interest (for figures & enrichement analysis)
	
	## Analysis Parameters
	Float? tol
	Float? var
	String? zscore_mode

	## Module Toggles
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
					label="${label}-so_nmf-${pair.left}",
					ome_pairs=[pair],
					gene_set_database=gene_set_database,
					yaml_file=yaml_file,
					groups_file=groups_file
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
				yaml_file=yaml_file,
				groups_file=groups_file
		}
	}

	# Sankey Diagrams
	if ( run_so_nmf && run_sankey ){ # if we ran so_nmf && want to run sankey
		call sankey_wdl.panoply_sankey_workflow as nmf_sankey {
			input:
				label = label,

				annot_files = so_nmf.nmf_membership,		# array of so-NMF results
				annot_file_labels = ome_labels,				# array of ome labels

				annot_file_primary = mo_nmf.nmf_membership, # single file with mo-NMF results
				annot_label_primary = "Multiomic",			# label for mo-NMF data

				annot_of_comparison="NMF.consensus",		# column for analysis
				annot_prefix=""								# no annot prefix to add
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
			sankey_report = nmf_sankey.sankey_report_file
	}


	output {
		File nmf_results=nmf_assemble.nmf_results
		File nmf_reports=nmf_assemble.nmf_reports
		# File? sankey_report=nmf_sankey.sankey_report_file
	}
}
