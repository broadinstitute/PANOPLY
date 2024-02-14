#
# Copyright (c) 2024 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_preprocess_gct/versions/1/plain-WDL/descriptor" as preprocess_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_ssgsea/versions/7/plain-WDL/descriptor" as ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_ssgsea_report/versions/8/plain-WDL/descriptor" as ssgsea_report_wdl

################################################
##  workflow: panoply_preprocess_gct + panoply_ssgsea + panoply_ssgsea_report
workflow panoply_ssgsea_workflow {

	File input_ds
	File gene_set_database
	File yaml_file
	String output_prefix

	Boolean preprocess_gct # toggle to convert GCT to gene-centric / single-site-centric


	## Preprocess GCT (optional) // Convert GCT to gene-centric or single-site-centric
	if (preprocess_gct) {
		call preprocess_wdl.panoply_preprocess_gct as preprocess {
		input:
			input_ds = input_ds,
			yaml_file = yaml_file
		}
	}

	## Run ssGSEA
	call ssgsea_wdl.panoply_ssgsea as ssgsea {
	input:
		input_ds = if defined(preprocess.result) then preprocess.result else input_ds,
		gene_set_database = gene_set_database,
		output_prefix = output_prefix,
		yaml_file = yaml_file
	}

	## Generate ssGSEA Report
    call ssgsea_report_wdl.panoply_ssgsea_report as ssgsea_report {
		input:
			tarball=ssgsea.results,
			cfg_yaml=yaml_file,
			label=output_prefix
	}

	output {
		File results=ssgsea.results
		File report=ssgsea_report.report
	}
}