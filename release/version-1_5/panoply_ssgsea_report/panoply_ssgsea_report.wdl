#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_ssgsea_report {

  File tarball
  File cfg_yaml
  Boolean? is_ptmsigdb
  String label

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-ssgsea.r -t ${tarball} -l ${label} -y ${cfg_yaml} -z /home/pgdac/src/ -p ${default=FALSE is_ptmsigdb}
  }

  output {
    File report = "report_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/panoply_ssgsea_report:1_5"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Karsten Krug"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_ssgsea_report_workflow {
  call panoply_ssgsea_report
}

