#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_harmonize/versions/2/plain-WDL/descriptor" as harmonize_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_sampleqc/versions/2/plain-WDL/descriptor" as sampleqc_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_sampleqc_report/versions/1/plain-WDL/descriptor" as sampleqc_report_wdl

workflow panoply_sampleqc_workflow {
  File inputData # gct if standalone='true' or tar if 'false'
  String type 
  File rnaExpr
  File cnaExp
  String job_identifier
  String standalone = "true"
  File yaml = '/prot/proteomics/Projects/PGDAC/src/master-parameters.yaml'

  Float? na_max
  String? duplicate_gene_policy
  String? gene_id_col

  call harmonize_wdl.panoply_harmonize {
    input:
      inputData = inputData,
      rnaExpr = rnaExpr,
      cnaExpr = cnaExpr,
      standalone = standalone,
      analysisDir = job_identifier,
      type = type,
      yaml = yaml,
      na_max = na_max,
      duplicate_gene_policy = duplicate_gene_policy,
      gene_id_col = gene_id_col
  }

  call sampleqc_wdl.panoply_sampleqc {
    input:
      tarball = panoply_harmonize.outputs,
      type = type,
      yaml = yaml
  }

  call sampleqc_report_wdl.panoply_sampleqc_report {
    input:
      tarball = panoply_sampleqc.outputs,
      type = type,
      label = job_identifier,
      tmpDir = "tmp"
  }

  output {
  	File harmonize_tar = panoply_harmonize.outputs
  	File sampleqc_tar = panoply_sampleqc.outputs
  	File sampleqc_report = panoply_sampleqc_report.report
  }
}