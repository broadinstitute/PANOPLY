#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_normalize_ms_data_report {
  File tarball
  String label
  String type
  String tmpDir
  File yaml

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    # Find the flag for normalize.proteomics in the yaml:
    cfg=${yaml}
    echo "library(yaml);yaml=read_yaml('$cfg');norm=yaml[['normalize.proteomics']];writeLines(as.character(norm), con='norm.txt')" > cmd.r
    Rscript cmd.r
    norm=`cat norm.txt`
    
    if [ $norm = FALSE ]; then
      echo 'no normalization performed' > "norm_"${label}".html"
    else
      Rscript /home/pgdac/src/rmd-normalize.r ${tarball} ${label} ${type} ${tmpDir}
    fi
  }

  output {
    File report = "norm_" + "${label}" + ".html"
  }

  runtime {
    docker : "broadcptacdev/panoply_normalize_ms_data_report:latest"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Karsten Krug"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_normalize_ms_data_report_workflow {
  call panoply_normalize_ms_data_report
}

