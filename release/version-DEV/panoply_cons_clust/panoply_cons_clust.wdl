#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cons_clust {
  File inputData   # output from panoply_harmonize or panoply_filter
  String type
  File? groupsFile
  File yaml
  String standalone
  String? analysisDir
  Int? clustering_sd_threshold
  Float? clustering_na_threshold
  String outFile = "panoply_cluster-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    codeDir="/prot/proteomics/Projects/PGDAC/src"
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module cons_clust \
    --master_yaml ${yaml} \
    ${"--clustering_sd_threshold " + clustering_sd_threshold} \
    ${"--clustering_na_threshold " + clustering_na_threshold} 
    echo ${type}
    
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh cluster \
            -i ${inputData} \
            -t ${type} \
            -c $codeDir \
            -o ${outFile} \
            -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
            ${"-g " + groupsFile} \
            -y "final_output_params.yaml";
   else
     /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh cluster \
            -f ${inputData} \ 
            -t ${type} \ 
            -c $codeDir \
            -r ${analysisDir} \
            -o ${outFile} \
            -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
            ${"-g " + groupsFile} \
            -y "final_output_params.yaml"
  fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptacdev/panoply_cons_clust:DEV"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu : select_first ([num_threads, 8]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani, Karsten Krug"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_cons_clust_workflow {
  call panoply_cons_clust
}

