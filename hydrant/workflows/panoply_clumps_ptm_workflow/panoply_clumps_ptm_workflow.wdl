#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm_diffexp/versions/4/plain-WDL/descriptor" as diffexp_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm_mapping/versions/12/plain-WDL/descriptor" as mapping_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm/versions/19/plain-WDL/descriptor" as analysis_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm_postprocess/versions/4/plain-WDL/descriptor" as postprocess_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm_report/versions/4/plain-WDL/descriptor" as report_wdl

################################################
##  workflow: panoply_clumps_ptm_diffexp + panoply_clumps_ptm_mapping + panoply_clumps_ptm
workflow panoply_clumps_ptm_workflow {
	# PTM GCT files; must include at least one
	File? pSTY_gct
	File? acK_gct
	File? ubK_gct

	String PDB_ref_bucket				# Google-Cloud Bucket with PDB Directory split into tarfiles

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
				PDB_ref_bucket = PDB_ref_bucket,
				yaml_file = yaml_file,
				output_prefix = output_prefix
		}
	}

  	scatter (diff_exp in diffexp.diff_exp_files) {
		call analysis_wdl.panoply_clumps_ptm as analysis {
		    input:
		    	diff_exp_file = diff_exp,
		    	var_sites_file = "${if defined(mapping_file) then mapping_file else mapping.filt_results}",
		        PDB_ref_bucket = PDB_ref_bucket,
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
	    	postprocess_results = postprocess.results, 			# array of results file from postprocess module
	    	label = output_prefix
	}

	output{
	    Array[File] clumps_ptm_results = analysis.results
	    File clumps_ptm_report = report.report
	}
}
