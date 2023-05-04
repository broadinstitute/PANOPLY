#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_so_nmf_sankey/versions/3/plain-WDL/descriptor" as sankey_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_so_nmf_sankey_report/versions/3/plain-WDL/descriptor" as so_nmf_report_wdl


workflow panoply_so_nmf_sankey_workflow {
  File so_nmf_tar
  File? mo_nmf_tar
  String label

  ## generate Sankey Diagrams comparing clustering results between -omes
  call sankey_wdl.panoply_so_nmf_sankey as nmf_sankey {
    input:
      so_nmf_tar = so_nmf_tar,
      mo_nmf_tar = mo_nmf_tar,
      label = label
  }

  ## generate report with sankey diagrams
  call so_nmf_report_wdl.panoply_so_nmf_sankey_report as nmf_sankey_report {
    input:
      sankey_tar = nmf_sankey.tar_out,
      label = label
  }
  


  output {
    File sankey_tar = nmf_sankey.tar_out
    File sankey_report = nmf_sankey_report.report_out
  }
 }