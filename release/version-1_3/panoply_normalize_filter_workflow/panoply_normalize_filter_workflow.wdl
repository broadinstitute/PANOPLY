#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_normalize_ms_data/versions/25/plain-WDL/descriptor" as normalize_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_filter/versions/2/plain-WDL/descriptor" as filter_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_normalize_ms_data_report/versions/17/plain-WDL/descriptor" as normalize_report_wdl

workflow panoply_normalize_filter_workflow {
	File input_pome
	String ome_type
	String job_identifier
	File yaml
	String? normalizeProteomics # "true" or "false"
	String? filterProteomics # "true" or "false"

	String? geneIdCol
  	String? proteinIdCol
  	String? proteinIdType
	String? combineReplicates
	Int? ndigits
  	Float? naMax
  	String? noNA
  	Float? sdFilterThreshold
  	

	call normalize_wdl.panoply_normalize_ms_data {
    	input:
      		inputData = input_pome,
      		type = ome_type,
      		standalone = "true",
      		analysisDir = job_identifier,
      		yaml = yaml,
      		ndigits = ndigits,
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

	call filter_wdl.panoply_filter {
    	input:
      		inputData = panoply_normalize_ms_data.outputs,
      		type = ome_type,
      		standalone = "true",
      		analysisDir = job_identifier,
      		yaml = yaml,
      		geneIdCol = geneIdCol,
      		proteinIdCol = proteinIdCol,
      		proteinIdType = proteinIdType,
      		filterProteomics=filterProteomics,
      		combineReplicates = combineReplicates,
      		ndigits = ndigits,
      		naMax = naMax,
      		noNA = noNA,
      		sdFilterThreshold = sdFilterThreshold
  	}

  	output {
  		File filtered_data_table = panoply_filter.outputs
  		File filtered_tar = panoply_filter.output_tar
  		File normalize_report = panoply_normalize_ms_data_report.report

  	}
  }