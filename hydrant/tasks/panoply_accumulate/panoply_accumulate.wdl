#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_accumulate {
  File input_tar
  String? output_tar = "panoply_contrasts.tar"
  String module
  String? analysisDir = "input_tarball"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/accumulate.sh \
      -i ${input_tar} \
      -o ${output_tar} \
      -r ${analysisDir} \
      -m ${module}
  }

  output {
    File outputs = "${output_tar}"
    Array[File] list_gct = glob( "${analysisDir}/${module}/contrasts/*.gct" )
  }

  runtime {
    docker      : "broadcptacdev/panoply_accumulate:latest"
    memory      : select_first ([memory, 16]) + "GB"
    disks       : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu         : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email  : "rkothadi@broadinstitute.org"
  }
}

workflow panoply_accumulate_workflow {
  call panoply_accumulate
}
