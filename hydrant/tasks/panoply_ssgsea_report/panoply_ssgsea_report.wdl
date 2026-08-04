#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

task panoply_ssgsea_report {

    File tarball
    File cfg_yaml
    String label

    # Heatmap Parameters
    Float? fdr
    Int? top_n
    Boolean? cluster_rows
    String? ser_meth
    Boolean? split_by_prefix

    File? geneset_groups_file

    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/ssgsea-renderRMD.R -t ${tarball} -l ${label} \
      ${"-g " + geneset_groups_file} \
      ${if defined(split_by_prefix) then "-s ${if select_first([split_by_prefix]) then 'TRUE' else 'FALSE'}" else ""} \
      ${"-f " + fdr} ${"-n " + top_n} \
      ${if defined(cluster_rows) then "-c ${if select_first([cluster_rows]) then 'TRUE' else 'FALSE'}" else ""} ${"-m " + ser_meth} \
      -y ${cfg_yaml} -z /home/pgdac/src/
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
