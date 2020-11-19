#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_parse_sm_table/versions/12/plain-WDL/descriptor" as parse_sm_table
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_normalize_ms_data/versions/16/plain-WDL/descriptor" as normalize
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_normalize_ms_data_report/versions/10/plain-WDL/descriptor" as normalize_report

workflow panoply_process_SM_table {

  ## inputs
  String job_identifier
  String ome_type
  File sample_annotation
  File input_ssv
  File yaml


  call parse_sm_table.panoply_parse_sm_table as parse {
    input:
      SMtable = input_ssv,
      exptDesign = sample_annotation,
      analysisDir = job_identifier,
      type = ome_type,
      yaml = yaml
  }

  call normalize.panoply_normalize_ms_data as norm {
    input:
      inputData = parse.outputs, 
      normalizeProteomics = "true",
      standalone = "false",
      type = ome_type,
      analysisDir = job_identifier,
      yaml = yaml
  }

  call normalize_report.panoply_normalize_ms_data_report as report {
    input:
      tarball = norm.output_tar,
      label = job_identifier,
      type = ome_type,
      tmpDir = "tmp",
      yaml = norm.output_yaml
  }
  
  output {
    File output_tar = norm.outputs
    File output_report = report.report
  }

}
