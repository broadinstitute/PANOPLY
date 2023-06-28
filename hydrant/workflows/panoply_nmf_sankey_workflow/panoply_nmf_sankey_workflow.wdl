#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_sankey/versions/1/plain-WDL/descriptor" as sankey_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_sankey_report/versions/1/plain-WDL/descriptor" as sankey_report_wdl


workflow panoply_nmf_sankey_workflow {
  Array[File] so_nmf_membership
  File? mo_nmf_membership
  String label

  ## generate Sankey Diagrams comparing clustering results between -omes
  call sankey_wdl.panoply_nmf_sankey as nmf_sankey {
    input:
      so_nmf_membership = so_nmf_membership,
      mo_nmf_membership = mo_nmf_membership,
      label = label
  }

  ## generate report with sankey diagrams
  call sankey_report_wdl.panoply_nmf_sankey_report as nmf_sankey_report {
    input:
      sankey_tar = nmf_sankey.tar_out,
      label = label
  }
  


  output {
    File sankey_tar = nmf_sankey.tar_out
    File sankey_report = nmf_sankey_report.report_out
  }
 }