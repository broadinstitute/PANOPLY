#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_association/versions/22/plain-WDL/descriptor" as assoc_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_accumulate/versions/20/plain-WDL/descriptor" as accum_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea/versions/29/plain-WDL/descriptor" as ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_association_report/versions/20/plain-WDL/descriptor" as 	assoc_report_wdl

################################################
##  workflow: panoply_association + panoply_accumulate + panoply_ssgsea + panoply_association_report
workflow panoply_association_workflow {

  	String job_identifier
  	String ome_type
	String standalone

	## inputs
	File inputData
	File yaml
	File association_groups

	Float? sample_na_max
	Float? nmiss_factor
	String? duplicate_gene_policy
	String? gene_id_col

	String geneset_db


	call assoc_wdl.panoply_association as assoc {
    input: 
    	inputData = inputData,
    	groupsFile = association_groups,
    	type = ome_type,
    	standalone = standalone,
    	analysisDir = job_identifier,
    	yaml = yaml,
    	sample_na_max=sample_na_max,
    	nmiss_factor=nmiss_factor,
    	duplicate_gene_policy=duplicate_gene_policy,
    	gene_id_col=gene_id_col
	}

	call accum_wdl.panoply_accumulate as accumulate_assoc {
	input:
		input_tar = assoc.outputs,
		module = "association"
	}

	Array[File] list_gct_assoc = accumulate_assoc.list_gct
	scatter (f in list_gct_assoc){
		call ssgsea_wdl.panoply_ssgsea as ssgsea_assoc {
		input:
			input_ds = "${f}",
			gene_set_database = geneset_db,
			output_prefix = job_identifier,
			level = "gc",
			yaml_file = yaml
		}
	}


	call assoc_report_wdl.panoply_association_report as assoc_report {
	input:
		ssgsea_assoc_tars = ssgsea_assoc.results,
	  	master_yaml = yaml,
  		label = "${job_identifier}-full-results",
  		type = ome_type
	}

	output {
		File outputs=assoc.outputs
		File contrasts=accumulate_assoc.outputs
		Array[File] ssgsea_assoc_tars=ssgsea_assoc.results
		File report=assoc_report.report_out
	}
}