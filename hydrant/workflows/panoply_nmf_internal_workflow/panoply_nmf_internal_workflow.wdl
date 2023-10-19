#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_balance_omes/versions/1/plain-WDL/descriptor" as panoply_nmf_balance_omes_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf/versions/1/plain-WDL/descriptor" as panoply_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_report/versions/1/plain-WDL/descriptor" as panoply_nmf_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_ssgsea/versions/1/plain-WDL/descriptor" as panoply_ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_ssgsea_report/versions/2/plain-WDL/descriptor" as panoply_ssgsea_report_wdl

################################################
##  workflow: nmf_balance_omes + nmf + nmf_report + ssgsea + ssgsea_report
workflow panoply_nmf_internal_workflow {
	
	String label
	Array[Pair[String,File]]+ ome_pairs

	File gene_set_database
	File? yaml_file
	File? groups_file

	# balance filter
	Boolean? balance_omes

	Float? tol
	Float? var
	String? zscore_mode


	# Toggle Balance Module -- run if we have multi-omic data && balancing is on
	if (length(omes) > 1 && balance_omes) {
		call panoply_nmf_balance_omes_wdl.panoply_nmf_balance_omes {
			input:
				label=label,
				ome_pairs=ome_pairs,
				tol=tol,
				var=var,
				zscore_mode=zscore_mode
		}
	}

	call panoply_nmf_wdl.panoply_nmf {
		input:
			ome_pairs="${if defined(panoply_nmf_balance_omes.ome_pairs)
						 then panoply_nmf_balance_omes.ome_pairs
						 else ome_pairs}",
			output_prefix=label,
			yaml_file=yaml_file
	}
	call panoply_nmf_wdl.panoply_nmf_postprocessing {
		input:
			nmf_results=panoply_nmf.results,
			nclust=panoply_nmf.nclust,
			expr_comb=panoply_nmf.gct_comb,
			expr_comb_nn=panoply_nmf.gct_comb_nn,
			output_prefix=label,
			groups_file=groups_file,
			yaml_file=yaml_file
	}

	call panoply_nmf_report_wdl.panoply_nmf_report {
        input:
			tarball=panoply_nmf.results,
			label=label
	}

	call panoply_ssgsea_wdl.panoply_ssgsea {
		input:
			input_ds=panoply_nmf.feature_matrix_w,
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
		File nmf_tar=panoply_nmf.results # full output tar
		File nmf_membership=panoply_nmf.membership # .txt with membership results
		Int  nmf_nclust=panoply_nmf.nclust # number of membershipers
		File nmf_report=panoply_nmf_report.report

		File nmf_ssgsea_tar=panoply_ssgsea.results
		File nmf_ssgsea_report=panoply_ssgsea_report.report
		File nmf_balance_filter=panoply_nmf_balance_omes.pdf
	}
}
