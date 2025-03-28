#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm_diffexp/versions/2/plain-WDL/descriptor" as diffexp_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm_mapping/versions/1/plain-WDL/descriptor" as mapping_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_clumps_ptm/versions/2/plain-WDL/descriptor" as analysis_wdl

################################################
##  workflow: panoply_clumps_ptm_diffexp + panoply_clumps_ptm_mapping + panoply_clumps_ptm
workflow panoply_clumps_ptm_workflow {
	# PTM GCT files; must include at least one
	File? pSTY_gct
	File? acK_gct
	File? ubK_gct

	File PDB_DIR				# PDB Directory tarfile

	String? output_prefix
	File yaml_file

	call diffexp_wdl.panoply_clumps_ptm_diffexp as diffexp {
	    input:
	        pSTY_gct = pSTY_gct,
	        acK_gct = acK_gct,
	        ubK_gct = ubK_gct,
	        yaml_file = yaml_file,
	        output_prefix = output_prefix
	}

	call mapping_wdl.panoply_clumps_ptm_mapping as mapping {
	    input:
	        pSTY_gct = pSTY_gct,
	        acK_gct = acK_gct,
	        ubK_gct = ubK_gct,
	        PDB_DIR = PDB_DIR,
	        yaml_file = yaml_file,
	        output_prefix = output_prefix
	}

  	scatter (diff_exp in diffexp.diff_exp_files) {
		call analysis_wdl.panoply_clumps_ptm as analysis {
		    input:
		    	diff_exp_file = diff_exp,
		    	var_sites_file = mapping.filt_results,
		        PDB_DIR = PDB_DIR,
	        	yaml_file = yaml_file,
		        output_prefix = output_prefix
		}
	}

	output{
	    Array[File] results = analysis.results
	}
}
