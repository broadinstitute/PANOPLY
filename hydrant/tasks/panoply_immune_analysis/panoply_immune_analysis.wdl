#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_immune_analysis {
  File inputData
  String type
  String standalone
  File yaml
  String? analysisDir
  File? groupsFile
  Float? fdr
  Int? heatmapWidth
  Int? heatmapHeight

  String outFile = "panoply_immune_analysis-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    codeDir="/prot/proteomics/Projects/PGDAC/src"
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r --module immune_analysis --master_yaml ${yaml} ${"--immune_enrichment_subgroups " + groupsFile} ${"--immune_enrichment_fdr " + fdr} ${"--immune_heatmap_width " + heatmapWidth} ${"--immune_heatmap_height " + heatmapHeight}
    
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh immune \
                  -i ${inputData} \
                  -t ${type} \
                  -c $codeDir \
                  -o ${outFile} \
                  ${"-g " + groupsFile} \
                  ${"-z " + fdr} \
                  -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
                  -y "final_output_params.yaml";
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh immune \
                  -rna ${inputData} \
                  -t ${type} \
                  -c $codeDir \
                  -r ${analysisDir} \
                  -o ${outFile} \
                  ${"-g " + groupsFile} \
                  ${"-z " + fdr} \
                  -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
                  -y "final_output_params.yaml"
    fi
  }

  output {
    File outputs = "${outFile}"
    File yaml_file = "final_output_params.yaml"
  }

  runtime {
    docker : "broadcptacdev/panoply_immune_analysis:latest"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_immune_analysis_workflow {
    call panoply_immune_analysis 
}
