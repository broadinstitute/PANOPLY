#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_metaboanalyst/versions/5/plain-WDL/descriptor" as metaboanalyst_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_metaboanalyst_report/versions/4/plain-WDL/descriptor" as metaboanalyst_report_wdl

################################################
##  workflow: panoply_metaboanalyst + panoply_metaboanalyst_report
workflow panoply_metaboanalyst_workflow {
	File meta_gct
	File? omic_gct
	String? ome_type

	String output_prefix
	File groups_file
	File yaml_file

	call metaboanalyst_wdl.panoply_metaboanalyst {
	    input:
	    	meta_gct = meta_gct,
	    	omic_gct = omic_gct,
	    	ome_type = ome_type,
	        output_prefix = output_prefix,
	        groups_file = groups_file,
	        yaml_file = yaml_file
	}

	call metaboanalyst_report_wdl.panoply_metaboanalyst_report {
	    input:
	        metaboanalyst_results = panoply_metaboanalyst.results,
	        label = output_prefix
	}

	output{
	    File metaboanalyst_tar = panoply_metaboanalyst.results
	    File metaboanalyst_report = panoply_metaboanalyst_report.report
	}
}
