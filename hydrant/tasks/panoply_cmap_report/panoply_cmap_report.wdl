#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cmap_report {

  File cmap_results
  File cmap_ssgsea_results
  String label

  # # Heatmap Parameters
  # Float? fdr
  # Int? top_n
  # Boolean? cluster_rows
  # String? ser_meth
  # Boolean? split_by_prefix

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/cmap-renderRMD.R \
      -c ${cmap_results} -s ${cmap_ssgsea_results} \
      -l ${label} -z /home/pgdac/src/
  }

  output {
    File report = label + "_CMAP_report.html"
  }

  runtime {
    docker : "broadcptacdev/panoply_cmap_report:latest"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "C.M. Williams"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_cmap_report_workflow {
  call panoply_cmap_report
}

