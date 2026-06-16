#
# Copyright (c) 2024 The Broad Institute, Inc. All rights reserved.
#
import "https://raw.githubusercontent.com/broadinstitute/PANOPLY/issue-githubWDL/hydrant/tasks/panoply_preprocess_gct/panoply_preprocess_gct.wdl" as preprocess_wdl
import "https://raw.githubusercontent.com/broadinstitute/PANOPLY/issue-githubWDL/hydrant/tasks/panoply_ssgsea/panoply_ssgsea.wdl" as ssgsea_wdl
import "https://raw.githubusercontent.com/broadinstitute/PANOPLY/issue-githubWDL/hydrant/tasks/panoply_ssgsea_report/panoply_ssgsea_report.wdl" as ssgsea_report_wdl

################################################
##  workflow: panoply_preprocess_gct + panoply_ssgsea + panoply_ssgsea_report
workflow panoply_ssgsea_workflow {

	File input_ds
	File gene_set_database
	File yaml_file
	String output_prefix

	Boolean preprocess_gct # toggle to convert GCT to gene-centric / single-site-centric

	## parameters to create gene-centric or single-site-centric 
    ## GCT files for ssGSEA / PTM-SEA
	String? level
 	String? id_type
	String? id_type_out
	String? acc_type
	String? seqwin_col
	String? gene_col
	Boolean? humanize_gene
	String? SGT_col
	Boolean? loc
	String? mode
	String? mod_res
	String? mod_type
	
	## ssGSEA / PTM-SEA parameters below	
	String? sample_norm_type
	String? correl_type
	String? statistic
	String? output_score_type
	Float? weight
	Int? min_overlap
	String? tolerate_min_overlap_err # boolean value: should the WDL tolerate "not-enough-overlap" errors?
	Int? nperm
	Boolean? global_fdr

	## Preprocess GCT (optional) // Convert GCT to gene-centric or single-site-centric
	if (preprocess_gct) {
		call preprocess_wdl.panoply_preprocess_gct as preprocess {
		input:
			input_ds = input_ds,
			gene_col = gene_col,
			yaml_file = yaml_file,

			level = level,
			id_type = id_type,
			id_type_out = id_type_out,
			acc_type = acc_type,
			seqwin_col = seqwin_col,
			gene_col = gene_col,
			humanize_gene = humanize_gene,
			SGT_col = SGT_col,
			loc = loc,
			mode = mode,
			mod_res = mod_res,
			mod_type = mod_type
		}
	}

	## Run ssGSEA
	call ssgsea_wdl.panoply_ssgsea as ssgsea {
	input:
		input_ds = if defined(preprocess.result) then preprocess.result else input_ds,
		gene_set_database = gene_set_database,
		output_prefix = output_prefix,
		yaml_file = yaml_file,
		
		sample_norm_type = sample_norm_type,
		correl_type = correl_type,
		statistic = statistic,
		output_score_type = output_score_type,
		weight = weight,
		min_overlap = min_overlap,
		tolerate_min_overlap_err = tolerate_min_overlap_err, # boolean value: should the WDL tolerate "not-enough-overlap" errors?
		nperm = nperm,
		global_fdr = global_fdr
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