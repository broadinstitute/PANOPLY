#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_association/versions/9/plain-WDL/descriptor" as assoc_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_accumulate/versions/8/plain-WDL/descriptor" as  accum_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_ssgsea/versions/6/plain-WDL/descriptor" as panoply_ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_association_report/versions/6/plain-WDL/descriptor" as 	assoc_report_wdl

################################################
##  workflow: mo_nmf_pre + mo_nmf + mo_nmf_report + ssgsea + ssgsea_report
workflow panoply_mo_nmf_gct_workflow {

  	String job_identifier
  	String ome_type

	## inputs
	File input_pome
	File yaml

	Float? sample_na_max
	Float? nmiss_factor
	String? duplicate_gene_policy
	String? gene_id_col

	String geneset_db


	call assoc_wdl.panoply_association {
    input: 
    	inputData = ome,
    	groupsFile = association_groups,
    	type = ome_type,
    	standalone = "true",
    	yaml = yaml_file,
    	sample_na_max=sample_na_max,
    	nmiss_factor=nmiss_factor,
    	duplicate_gene_policy=duplicate_gene_policy,
    	gene_id_col=gene_id_col
	}

	call accum_wdl.panoply_accumulate as accumulate_assoc {
	input:
		input_tar = panoply_association.outputs,
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




	call assoc_report_wdl.panoply_association_report {
	input:
		input_tar = panoply_download.full,
	  	master_yaml = yaml,
  		label = "${job_identifier}-full-results",
  		type = ome_type
	}

	output {
		File nmf_clust=panoply_mo_nmf.results
		File nmf_clust_report=panoply_mo_nmf_report.report
		File nmf_ssgsea=panoply_ssgsea.results
		File nmf_ssgsea_report=panoply_ssgsea_report.report
		File nmf_balance_filter=panoply_mo_nmf_pre.pdf
	}
}
