#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm_diffexp/versions/2/plain-WDL/descriptor" as diffexp_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm_mapping/versions/9/plain-WDL/descriptor" as mapping_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm/versions/7/plain-WDL/descriptor" as analysis_wdl
#import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm_report/versions/1/plain-WDL/descriptor" as report_wdl

################################################
##  workflow: panoply_clumps_ptm_diffexp + panoply_clumps_ptm_mapping + panoply_clumps_ptm
workflow panoply_clumps_ptm_workflow {
	# PTM GCT files; must include at least one
	File? pSTY_gct
	File? acK_gct
	File? ubK_gct

	Array[File]+ PDB_DIR				# PDB Directory in (multiple) tarfiles

	File? mapping_file					# pre-generated mapping file, to skip mapping module

	String? output_prefix
	File yaml_file

	call diffexp_wdl.panoply_clumps_ptm_diffexp as diffexp {
		input:
			pSTY_gct = pSTY_gct,
			acK_gct = acK_gct,
			ubK_gct = ubK_gct,
			yaml_file = yaml_file,
			output_prefix = output_prefix
	}

    if ( !defined(mapping_file) ){ # only run if a mapping file is not provided
		call mapping_wdl.panoply_clumps_ptm_mapping as mapping {
			input:
				pSTY_gct = pSTY_gct,
				acK_gct = acK_gct,
				ubK_gct = ubK_gct,
				PDB_DIR = PDB_DIR,
				yaml_file = yaml_file,
				output_prefix = output_prefix
		}
	}

  	scatter (diff_exp in diffexp.diff_exp_files) {
		call analysis_wdl.panoply_clumps_ptm as analysis {
		    input:
		    	diff_exp_file = diff_exp,
		    	var_sites_file = "${if defined(mapping_file) then mapping_file else mapping.filt_results}",
		        PDB_DIR = PDB_DIR,
	        	yaml_file = yaml_file,
		        output_prefix = sub(basename(diff_exp), "_diff_exp\\.tsv$", "") # use diffexp file as output prefix (to keep annotation)
		}
	}

	#call report_wdl.panoply_clumps_ptm_report as report {
	#    input:
	#    	results = analysis.results 			# array of results file from analysis module
	#    	output_prefix = output_prefix
	#}

	output{
	    Array[File] results = analysis.results
	    #File report = report.report
	}
}
