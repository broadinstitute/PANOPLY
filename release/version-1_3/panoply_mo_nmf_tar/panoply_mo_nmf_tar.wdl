#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_mo_nmf/versions/18/plain-WDL/descriptor" as panoply_mo_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea/versions/27/plain-WDL/descriptor" as panoply_ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea_report/versions/19/plain-WDL/descriptor" as panoply_ssgsea_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_mo_nmf_report/versions/11/plain-WDL/descriptor" as panoply_mo_nmf_report_wdl

################################################
##  workflow: mo_nmf + mo_nmf_report + ssgsea + ssgsea_report
workflow panoply_mo_nmf_tar_workflow {
	
	File tar_file
	File gene_set_database
        File yaml_file
	String label

	call panoply_mo_nmf_wdl.panoply_mo_nmf {
		input:
			tar_file=tar_file,
			yaml_file=yaml_file
	}

	call panoply_mo_nmf_report_wdl.panoply_mo_nmf_report {
 	        input:
			tarball=panoply_mo_nmf.results,
			label=label
	}
	
	call panoply_ssgsea_wdl.panoply_ssgsea {
		input:
			input_ds=panoply_mo_nmf.feature_matrix_w,
			gene_set_database=gene_set_database,
			yaml_file=yaml_file,
			output_prefix=label,
 			mode="abs.max",
			weight=1,
			
	}

    	call panoply_ssgsea_report_wdl.panoply_ssgsea_report {
	        input:
		        tarball=panoply_ssgsea.results,
			cfg_yaml=yaml_file,
			label=label
		
	}

	output {
		File nmf_clust=panoply_mo_nmf.results
		File nmf_clust_report=panoply_mo_nmf_report.report
		File nmf_ssgsea=panoply_ssgsea.results
		File nmf_ssgsea_report=panoply_ssgsea_report.report
		}
}
