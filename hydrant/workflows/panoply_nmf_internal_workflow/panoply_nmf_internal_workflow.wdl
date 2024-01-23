#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_balance_omes/versions/7/plain-WDL/descriptor" as panoply_nmf_balance_omes_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf/versions/6/plain-WDL/descriptor" as panoply_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_postprocess/versions/5/plain-WDL/descriptor" as panoply_nmf_postprocess_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_report/versions/3/plain-WDL/descriptor" as panoply_nmf_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_ssgsea/versions/6/plain-WDL/descriptor" as panoply_ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_ssgsea_report/versions/8/plain-WDL/descriptor" as panoply_ssgsea_report_wdl

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

    # separate paired-array into array of labels and array of GCTs
    scatter (pairs in ome_pairs) {
        String ome_labels = pairs.left
        File ome_gcts = pairs.right
    }

	# Toggle Balance Module -- run if we have multi-omic data && balancing is on
	if (length(ome_pairs) > 1 && select_first([balance_omes, false])) { # false by default, if balance_omes not provided
		call panoply_nmf_balance_omes_wdl.panoply_nmf_balance_omes as balance {
			input:
				label=label,
				ome_gcts=ome_gcts,
				ome_labels=ome_labels,
				tol=tol,
				var=var,
				zscore_mode=zscore_mode
		}
	}

	call panoply_nmf_wdl.panoply_nmf as nmf {
		input:
			ome_gcts=if defined(balance.ome_gcts_balanced) then balance.ome_gcts_balanced else ome_gcts,
			ome_labels=ome_labels,
			output_prefix=label,
			yaml_file=yaml_file
	}
	call panoply_nmf_postprocess_wdl.panoply_nmf_postprocess as postprocess {
		input:
			nmf_results=nmf.results,
			nclust=nmf.nclust,
			output_prefix=label,
			groups_file=groups_file,
			yaml_file=yaml_file
	}

	call panoply_nmf_report_wdl.panoply_nmf_report {
        input:
			nmf_results=nmf.results,
			nclust=nmf.nclust,
            postprocess_tarball=postprocess.results,
			label=label
	}

	call panoply_ssgsea_wdl.panoply_ssgsea {
		input:
			input_ds=postprocess.feature_matrix_w,
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
		File nmf_results=nmf.results				## tar file output of NMF (combined expression input-GCTs & res.rank Rdata object)
		Int  nmf_nclust=nmf.nclust					## number of clusters
		File nmf_figures=postprocess.results		## tar file with figures & analysis & parameters
		File nmf_membership=postprocess.membership	## .tsv with membership results
		File nmf_report=panoply_nmf_report.report 	## report file

		File nmf_ssgsea_tar=panoply_ssgsea.results
		File nmf_ssgsea_report=panoply_ssgsea_report.report
        
		File? nmf_balance_filter=balance.pdf
		File? nmf_preprocess_figures=nmf.preprocess_figs
	}
}
