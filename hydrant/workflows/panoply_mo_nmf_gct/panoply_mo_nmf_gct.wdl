import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_mo_nmf_pre/versions/1/plain-WDL/descriptor" as panoply_mo_nmf_pre_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_mo_nmf/versions/2/plain-WDL/descriptor" as panoply_mo_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea/versions/9/plain-WDL/descriptor" as panoply_ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea_report/versions/1/plain-WDL/descriptor" as panoply_ssgsea_report_wdl

################################################
##  workflow: mo-nmf plus ssgsea
workflow panoply_mo_nmf_gct_workflow {

	File gene_set_database
        File yaml_file
	String label

	call panoply_mo_nmf_pre_wdl.panoply_mo_nmf_pre {
		input:
			label=label
	}

	call panoply_mo_nmf_wdl.panoply_mo_nmf {
		input:
			tar_file=panoply_mo_nmf_pre.tar,
			yaml_file=yaml_file
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
		File nmf_ssgsea=panoply_ssgsea.results
		File nmf_ssgsea_report=panoply_ssgsea_report.report
		}
}
