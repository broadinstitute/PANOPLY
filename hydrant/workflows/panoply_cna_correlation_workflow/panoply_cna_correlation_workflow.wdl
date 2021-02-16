#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_harmonize/versions/2/plain-WDL/descriptor" as harmonize_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_cna_setup/versions/2/plain-WDL/descriptor" as cna_setup_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_cna_correlation/versions/2/plain-WDL/descriptor" as cna_corr_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_cna_correlation_report/versions/1/plain-WDL/descriptor" as cna_corr_report_wdl

workflow panoply_cna_correlation_workflow {
  File inputData # gct if standalone='true' or tar if 'false'
  String type 
  File rnaExpr
  File cnaExp
  String job_identifier
  File? groupsFile
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

  call cna_setup_wdl.panoply_cna_setup {
    input:
      tarball = panoply_harmonize.outputs,
      groupsFile = groupsFile,
      type = type,
      yaml = yaml
  }

  call cna_corr_wdl.panoply_cna_correlation {
    input:
      tarball = panoply_cna_setup.outputs,
      type = type,
      yaml = yaml
  }

  call cna_corr_report_wdl.panoply_cna_correlation_report {
    input:
      tarball = panoply_cna_correlation.outputs,
      config_yaml = yaml,
      type = type,
      label = job_identifier,
      tmpDir = "tmp"
  }

  output {
  	File harmonize_tar = panoply_harmonize.outputs
  	File cna_correlation_tar = panoply_cna_correlation.outputs
  	File cna_corr_report = panoply_cna_correlation_report.report
  }
}


