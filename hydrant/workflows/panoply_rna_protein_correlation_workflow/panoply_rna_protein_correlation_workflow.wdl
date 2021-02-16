#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_rna_protein_correlation/versions/2/plain-WDL/descriptor" as rna_prot_corr_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_rna_protein_correlation_report/versions/1/plain-WDL/descriptor" as rna_corr_report_wdl

workflow panoply_rna_protein_correlation_workflow {
  File inputData
  File rnaExpr 
  String type
  String standalone = "true"
  String job_identifier 
  File yaml = '/prot/proteomics/Projects/PGDAC/src/master-parameters.yaml'

  call rna_prot_corr_wdl.panoply_rna_protein_correlation {
    input:
      inputData = inputData,
      type = type,
      rnaExpr = rnaExpr,
      analysisDir = job_identifier,
      standalone = standalone,
      yaml = yaml
  }

  call rna_corr_report_wdl.panoply_rna_protein_correlation_report {
    input:
      tarball = panoply_rna_protein_correlation.outputs,
      config_yaml = yaml,
      label = job_identifier,
      type = type,
      tmpDir = "tmp"
  }

  output {
  	File rna_prot_corr_output = panoply_rna_protein_correlation.outputs
  	File rna_corr_report = panoply_rna_protein_correlation_report.report
  }
}