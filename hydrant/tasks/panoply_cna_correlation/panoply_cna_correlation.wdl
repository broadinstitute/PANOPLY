#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cna_correlation {
  File tarball   # output from panoply_cna_setup
  String type
  Float? fdr_cna_corr
  File yaml

  String outFile = "panoply_cna_correlation-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    codeDir="/prot/proteomics/Projects/PGDAC/src"
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module cna_analysis \
    --master_yaml ${yaml} \
    ${"--fdr_cna_corr " + fdr_cna_corr}
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CNAcorr -i ${tarball} -t ${type} -c $codeDir -o ${outFile} -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" -y "final_output_params.yaml"
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptacdev/panoply_cna_setup:latest"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "proteogenomics@broadinstitute.org"
  }
}


workflow panoply_cna_correlation_workflow {
  call panoply_cna_correlation
}

