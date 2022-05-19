#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_normalize_ms_data/versions/20/plain-WDL/descriptor" as normalize_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_normalize_ms_data_report/versions/14/plain-WDL/descriptor" as normalize_report_wdl

workflow panoply_normalize_ms_data_workflow {
	File input_pome
	String ome_type
	String job_identifier
	File yaml
	String? normalizeProteomics # "true" or "false"

	Int? ndigits
  	Float? na_max
  	String? gene_id_col
  	Float? sd_filter_threshold
  	Float? min_numratio_fraction
  	

	call normalize_wdl.panoply_normalize_ms_data {
    	input:
      		inputData = input_pome,
      		type = ome_type,
      		standalone = "true",
      		analysisDir = job_identifier,
      		yaml = yaml,
      		ndigits = ndigits,
      		naMax = na_max,
      		geneIdCol=gene_id_col,
      		sdFilterThreshold=sd_filter_threshold,
      		minNumratioFraction=min_numratio_fraction,
      		normalizeProteomics=normalizeProteomics
  	}

  	call normalize_report_wdl.panoply_normalize_ms_data_report {
    	input:
      		tarball = panoply_normalize_ms_data.output_tar,
      		label = job_identifier,
      		type = ome_type,
      		tmpDir = "tmp",
      		yaml = panoply_normalize_ms_data.output_yaml
  	}

  	output {
  		File normalized_data_table = panoply_normalize_ms_data.outputs
  		File normalized_tar = panoply_normalize_ms_data.output_tar
  		File normalize_report = panoply_normalize_ms_data_report.report

  	}
  }