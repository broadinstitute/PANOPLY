#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cna_correlation_report {
  File tarball
  File config_yaml
  
  String label
  String type
  String tmpDir

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    cfg=${config_yaml}
    echo "library(yaml);yaml=read_yaml('$cfg');fdr=yaml[['panoply_cna_correlation_report']]\$fdr;writeLines(as.character(fdr), con='fdr.txt')" > cmd.r
    Rscript cmd.r
    fdr=`cat fdr.txt`
    Rscript /home/pgdac/src/rmd-cna-analysis.r ${tarball} ${label} ${type} $fdr ${tmpDir}
  }

  output {
    File report = "cna-analysis_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/panoply_cna_correlation_report:1_2"
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

workflow panoply_cna_correlation_report_workflow {
  call panoply_cna_correlation_report
}
