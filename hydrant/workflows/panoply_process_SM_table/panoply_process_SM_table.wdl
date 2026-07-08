#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
version 1.0

import "../../tasks/panoply_parse_sm_table/panoply_parse_sm_table.wdl" as parse_sm_table
import "../../tasks/panoply_normalize_ms_data/panoply_normalize_ms_data.wdl" as normalize
import "../../tasks/panoply_normalize_ms_data_report/panoply_normalize_ms_data_report.wdl" as normalize_report

workflow panoply_process_SM_table {

  ## inputs
  input {
    String job_identifier
    String ome_type
    File sample_annotation
    File input_ssv
    File yaml
  }

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
