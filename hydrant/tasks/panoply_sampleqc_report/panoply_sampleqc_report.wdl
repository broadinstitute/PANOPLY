#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_sampleqc_report {
  File tarball
  String label
  String type
  String tmpDir

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-sample-qc.r ${tarball} ${label} ${type} ${tmpDir}
  }

  output {
    File report = "sample-qc_" + label + ".html"
  }

  runtime {
    docker : "broadcptacdev/panoply_sampleqc_report:latest"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_sampleqc_report_workflow {
	call panoply_sampleqc_report
}
