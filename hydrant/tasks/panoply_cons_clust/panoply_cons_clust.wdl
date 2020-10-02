#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cons_clust {
  File tarball   # output from pgdac_harmonize or pgdac_normalize_ms_data
  String type
  File? groupsFile
  String? subType
  File yaml
  Int? clustering_sd_threshold
  Float? clustering_na_threshold
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_cluster-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module cons_clust \
    --master_yaml ${yaml} \
    ${"--clustering_sd_threshold " + clustering_sd_threshold} \
    ${"--clustering_na_threshold " + clustering_na_threshold} 
    echo ${type}
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh cluster -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-m " + subType} -p "config-custom.r" ${"-g " + groupsFile}
   
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptacdev/panoply_cons_clust:latest"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu : select_first ([num_threads, 8]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani, Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}

workflow panoply_cons_clust_workflow {
  call panoply_cons_clust
}

