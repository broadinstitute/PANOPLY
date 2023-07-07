#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_sankey/versions/1/plain-WDL/descriptor" as sankey_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_sankey_report/versions/1/plain-WDL/descriptor" as sankey_report_wdl


workflow panoply_sankey_workflow {
  Array[File]+ annot_files          # annotation file(s)
  Array[String]+ annot_file_labels  # datatype(s) / label(s) for the provided annotation file(s)

  File? annot_file_primary          # annotation file which should be centered / highlighted in sankey-diagrams
  String? annot_label_primary       # corresponding label

  String annot_of_comparison        # annotation column which should be used for sankey comparisons (e.g. "NMF.consensus")
  String annot_prefix               # prefix to prepent to annotation values (e.g. "C" -> C1 C2 C2, instead of 1 2 3)


  String label

  ## generate Sankey Diagrams comparing clustering results between -omes
  call sankey_wdl.panoply_sankey as sankey {
    input:
      annot_files=annot_files,
      annot_file_labels=annot_file_labels,

      annot_file_primary=annot_file_primary,
      annot_label_primary=annot_label_primary,

      annot_column=annot_of_comparison,
      annot_prefix=annot_prefix,

      label = label
  }

  ## generate report with sankey diagrams
  call sankey_report_wdl.panoply_sankey_report as sankey_report {
    input:
      sankey_tar=sankey.tar_out,
      annot_of_comparison=annot_of_comparison,
      primary_dataype_label=annot_label_primary,
      
      label=label
  }
  


  output {
    File sankey_tar = sankey.tar_out
    File sankey_report_file = sankey_report.report_out
  }
 }