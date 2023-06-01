#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_so_nmf_pre/versions/2/plain-WDL/descriptor" as panoply_so_nmf_pre_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_mo_nmf/versions/18/plain-WDL/descriptor" as panoply_mo_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_so_nmf_assemble_drivers/versions/3/plain-WDL/descriptor" as panoply_so_nmf_assemble_drivers
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_mo_nmf_report/versions/11/plain-WDL/descriptor" as panoply_mo_nmf_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea/versions/26/plain-WDL/descriptor" as panoply_ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea_report/versions/16/plain-WDL/descriptor" as panoply_ssgsea_report_wdl

################################################
##  workflow: so_nmf_pre + mo_nmf + mo_nmf_report + ssgsea + ssgsea_report
workflow panoply_so_nmf_gct_workflow {

  File gene_set_database
  File yaml_file
  String label
    
  File ome
  String ome_type

	call panoply_so_nmf_pre_wdl.panoply_so_nmf_pre {
		input:
			label=label,
            ome=ome,
            ome_type=ome_type
	}

	call panoply_mo_nmf_wdl.panoply_mo_nmf {
		input:
			tar_file=panoply_so_nmf_pre.tar,
			yaml_file=yaml_file,
            output_prefix="results_nmf-${label}"
	}
    
    
	call panoply_so_nmf_assemble_drivers.panoply_so_nmf_assemble_drivers {
		input:
			nmf_tar=panoply_mo_nmf.results,
            ome_type=ome_type
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
			output_prefix="results_ssgsea-${label}",
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
        Pair[String, File] nmf_drivers=panoply_so_nmf_assemble_drivers.driver_pair
	}
}
