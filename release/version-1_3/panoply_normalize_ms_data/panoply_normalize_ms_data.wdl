#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_normalize_ms_data {
  File inputData
  String type
  String standalone
  String analysisDir
  File yaml
  String? normalizeProteomics
  String? normMethod
  String? altMethod
  Int? ndigits

  String outTar = "panoply_normalize_ms_data-output.tar"
  String outTable = "normalized_table-output.gct"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command  <<<
    set -euo pipefail
    
    codeDir="/prot/proteomics/Projects/PGDAC/src"
    dataDir="/prot/proteomics/Projects/PGDAC/data"
    
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
      --module normalize_ms_data \
      --master_yaml ${yaml} \
      ${"--normalize_proteomics " + normalizeProteomics} \
      ${"--norm_method " + normMethod} \
      ${"--alt_method " + altMethod} \
      ${"--ndigits " + ndigits} \
    
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh normalize \
            -i ${inputData} \
            -t ${type} \
            -c $codeDir \
            -o ${outTar} \
            -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
            -y "final_output_params.yaml"
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh normalize \
            -a ${inputData} \
            -r ${analysisDir} \
            -t ${type} \
            -c $codeDir \
            -d $dataDir \
            -o ${outTar} \
            -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
            -y "final_output_params.yaml"
    fi
    # Grab the normalized gct to set as output with appropriate name
    outGCT=`find ${analysisDir}/normalized-data -type f -iname "*-ratio-norm.gct"`
    outTableName=${type}-${outTable} 
    cp $outGCT $outTableName
    
  >>>

  output {
    File outputs = "${type}-${outTable}"
    File output_tar = "${outTar}"
    File output_yaml = "final_output_params.yaml"
  }

  runtime {
    docker : "broadcptac/panoply_normalize_ms_data:1_3"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_normalize_ms_data_workflow {

  call panoply_normalize_ms_data

}