#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
version 1.1

task panoply_ssgsea_report {
  input {

    File tarball
    File cfg_yaml
    String label

    # Heatmap Parameters
    # Defaults are read from the 'panoply_ssgsea_report' section of cfg_yaml;
    # any parameter set here overrides the yaml value.
    Float? fdr
    Int? top_n
    Boolean? cluster_rows
    Boolean? split_by_prefix

    File? geneset_groups_file

    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions
  }

  command {
    set -euo pipefail
    # (parameters that are not set are omitted, so ssgsea-renderRMD.R takes them from cfg_yaml)
    Rscript /home/pgdac/src/ssgsea-renderRMD.R \
      --tarball ${tarball} --label ${label} --yaml ${cfg_yaml} \
      ${"--fdr " + fdr} \
      ${"--top_n " + top_n} \
      ${"--cluster_rows " + cluster_rows} \
      ${"--split_by_prefix " + split_by_prefix} \
      ${"--geneset_groups_file " + geneset_groups_file}
  }

  output {
    File report = label + "_ssGSEA_rmd.html"
  }

  runtime {
    docker : "broadcptacdev/panoply_ssgsea_report:latest"
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

