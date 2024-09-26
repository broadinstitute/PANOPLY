#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_filter {
  File inputData
  String type
  String standalone
  String analysisDir
  File yaml
  String? geneIdCol
  String? proteinIdCol
  String? proteinIdType
  String? filterProteomics
  String? separateQCTypes
  String? combineReplicates
  Int? ndigits
  Float? naMax
  String? noNA
  Float? sdFilterThreshold

  String outTar = "panoply_filter-output.tar"
  String outTable = "filtered_table-output.gct"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command  <<<
    set -euo pipefail
    
    codeDir="/prot/proteomics/Projects/PGDAC/src"
    dataDir="/prot/proteomics/Projects/PGDAC/data"
    
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
      --module filter \
      --master_yaml ${yaml} \
      ${"--filter_proteomics " + filterProteomics} \
      ${"--gene_id_col " + geneIdCol} \
      ${"--protein_id_col " + proteinIdCol} \
      ${"--protein_id_type " + proteinIdType} \
      ${"--combine_replicates " + combineReplicates} \
      ${"--separate_qc_types " + separateQCTypes} \
      ${"--ndigits " + ndigits} \
      ${"--na_max " + naMax} \
      ${"--no_na " + noNA} \
      ${"--sd_filter_threshold " + sdFilterThreshold}
    

    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh filter \
            -i ${inputData} \
            -t ${type} \
            -c $codeDir \
            -o ${outTar} \
            -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
            -y "final_output_params.yaml"
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh filter \
            -n ${inputData} \
            -r ${analysisDir} \
            -t ${type} \
            -c $codeDir \
            -d $dataDir \
            -o ${outTar} \
            -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
            -y "final_output_params.yaml"
    fi

    # Grab the filtered gct to set as output with appropriate name
    outGCT=`find ${analysisDir}/filtered-data -type f -iname "*-ratio-norm-filt.gct"` # grab filtered file
    outTableName=${type}-${outTable} 
    cp $outGCT $outTableName
  >>>

  output {
    File outputs = "${type}-${outTable}"
    File output_tar = "${outTar}"
    File output_yaml = "final_output_params.yaml"
  }

  runtime {
    docker : "broadcptac/panoply_filter:1_5"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "C.M. Williams"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_filter_workflow {

  call panoply_filter
}