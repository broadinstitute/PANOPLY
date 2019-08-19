task pgdac_normalize_ms_data {
  File inputData
  String type
  String standalone
  String? analysisDir
  String? subType
  File? params

  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_normalize_ms_data-output.tar"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh normalize \
              -i ${inputData} \
              -t ${type} \
              -c ${codeDir} \
              -o ${outFile} \
              ${"-m " + subType} \
              ${"-p " + params};
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh normalize \
              -a ${inputData} \
              -r ${analysisDir} \
              -t ${type} \
              -c ${codeDir} \
              -d ${dataDir} \
              -o ${outFile} \
              ${"-m " + subType} \
              ${"-p " + params};
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_normalize_ms_data:1"
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

workflow pgdac_normalize_ms_data_workflow {
  File inputData
  String dataType
  String standalone
  String? analysisDir

  call pgdac_normalize_ms_data {
    input:
      inputData=inputData,
      type=dataType,
      standalone=standalone,
      analysisDir=analysisDir
  }
}
