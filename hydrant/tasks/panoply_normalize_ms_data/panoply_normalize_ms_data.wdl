task panoply_normalize_ms_data {
  File inputData
  String type
  String standalone
  String? analysisDir
  String? subType
  File yaml
  String? normalizeProteomics
  String? normMethod
  String? altMethod
  Int? ndigits
  Float? na_max
  String? gene_id_col
  Float? sd_filter_threshold
  Float? min_numratio_fraction

  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outTar = "panoply_normalize_ms_data-output.tar"
  String outTable = "normalized_table-output.gct"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command  <<<
    set -euo pipefail
    
    if [ ${normalizeProteomics} ]; then
      if [ ${normalizeProteomics} = "FALSE" ]; then
        norm=FALSE
      fi
      if [ ${normalizeProteomics} = "TRUE" ]; then
        norm=TRUE
      fi
    else
      # Find the flag for normalize.proteomics in the yaml:
      cfg=${yaml}
      echo "library(yaml);yaml=read_yaml('$cfg');norm=yaml[['normalize.proteomics']];writeLines(as.character(norm), con='norm.txt')" > cmd.r
      Rscript cmd.r
      norm=`cat norm.txt`
    fi
    
    # If not normalizing input is the normalized tab and becomes the output, else run normalizing:
    if [ $norm = FALSE ]; then
      cp ${inputData} ${type}-${outTable}
      tar -c -f ${outTar} ${type}-${outTable}
    else
      Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
      --module normalize_ms_data \
      --master_yaml ${yaml} \
      ${"--norm_method " + normMethod} \
      ${"--alt_method " + altMethod} \
      ${"--ndigits " + ndigits} \
      ${"--gene_id_col " + gene_id_col} \
      ${"--na_max " + na_max} \
      ${"--sd_filter_threshold " + sd_filter_threshold} \
      ${"--min_numratio_fraction " + min_numratio_fraction}
      
      if [[ ${standalone} = false ]]; then
        /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh normalize \
              -i ${inputData} \
              -t ${type} \
              -c ${codeDir} \
              -o ${outTar} \
              ${"-m " + subType} \
              -p "config-custom.r"
      else
        /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh normalize \
              -a ${inputData} \
              -r ${analysisDir} \
              -t ${type} \
              -c ${codeDir} \
              -d ${dataDir} \
              -o ${outTar} \
              ${"-m " + subType} \
              -p "config-custom.r"
      fi
      # Grab the norm/filtered gct to set as output with appropriate name
      outGCT=`find ${analysisDir}/normalized-data -type f -iname "*-ratio-norm-NArm*"`
      outTableName=${type}-${outTable} 
      cp $outGCT $outTableName
    fi
  >>>

  output {
    File outputs = "${type}-${outTable}"
    File output_tar = "${outTar}"
  }

  runtime {
    docker : "broadcptacdev/panoply_normalize_ms_data:latest"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email : "rkothadi@broadinstitute.org"
  }
}

workflow panoply_normalize_ms_data_workflow {
  File inputData
  String dataType
  String standalone
  String? analysisDir
  File yaml
  Int? ndigits
  String? gene_id_col
  Float? na_max
  Float? sd_filter_threshold
  Float? min_numratio_fraction

  call panoply_normalize_ms_data {
    input:
      inputData=inputData,
      type=dataType,
      standalone=standalone,
      analysisDir=analysisDir,
      yaml=yaml,
      ndigits=ndigits,
      gene_id_col=gene_id_col,
      na_max=na_max,
      sd_filter_threshold=sd_filter_threshold,
      min_numratio_fraction=min_numratio_fraction
  }
}