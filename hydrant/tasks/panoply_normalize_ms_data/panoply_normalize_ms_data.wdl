task panoply_normalize_ms_data {
  File inputData
  String type
  String standalone
  String? analysisDir
  String? subType
  File yaml
  String? normMethod
  String? altMethod
  Int? ndigits
  Float? naMax
  Float? sampleNaMax
  Float? minNumratioFraction
  Float? nmissFactor
  Float? sdFilterThreshold
  String? duplicateGenePolicy

  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_normalize_ms_data-output.tar"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module normalize_ms_data \
    --master_yaml ${yaml} \
    ${"--norm_method " + normMethod} \
    ${"--alt_method " + altMethod} \
    ${"--ndigits " + ndigits} \
    ${"--na_max " + naMax} \
    ${"--sample_na_max " + sampleNaMax} \
    ${"--min_numratio_fraction " + minNumratioFraction} \
    ${"--nmiss_factor " + nmissFactor} \
    ${"--sd_filter_threshold " + sdFilterThreshold} \
    ${"--duplicate_gene_policy " + duplicateGenePolicy}
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh normalize \
              -i ${inputData} \
              -t ${type} \
              -c ${codeDir} \
              -o ${outFile} \
              ${"-m " + subType} \
              -p "config-custom.r";
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh normalize \
              -a ${inputData} \
              -r ${analysisDir} \
              -t ${type} \
              -c ${codeDir} \
              -d ${dataDir} \
              -o ${outFile} \
              ${"-m " + subType} \
              -p "config-custom.r";
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/panoply_normalize_ms_data:dev"
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
  String? normMethod
  String? altMethod
  Int? ndigits
  Float? naMax
  Float? sampleNaMax
  Float? minNumratioFraction
  Float? nmissFactor
  Float? sdFilterThreshold
  String? duplicateGenePolicy

  call panoply_normalize_ms_data {
    input:
      inputData=inputData,
      type=dataType,
      standalone=standalone,
      analysisDir=analysisDir,
      yaml=yaml,
      normMethod=normMethod,
      altMethod=altMethod,
      ndigits=ndigits,
      naMax=naMax,
      sampleNaMax=sampleNaMax,
      minNumratioFraction=minNumratioFraction,
      nmissFactor=nmissFactor,
      sdFilterThreshold=sdFilterThreshold,
      duplicateGenePolicy=duplicateGenePolicy
  }
}
