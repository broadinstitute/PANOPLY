task pgdac_harmonize {
  File inputData
  File rnaExpr
  File cnaExpr
  String type
  String connected
  String? analysisDir
  String? subType
  File? params

  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"
  String outFile = "pgdac_harmonize-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    if [[ ${connected} = true ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize \
                  -i ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -d ${dataDir} \
                  -rna ${rnaExpr} \
                  -cna ${cnaExpr} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize \
                  -f ${inputData} \
                  -r ${analysisDir} \
                  -t ${type} \
                  -c ${codeDir} \
                  -d ${dataDir} \
                  -rna ${rnaExpr} \
                  -cna ${cnaExpr} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_harmonize:1"
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

workflow pgdac_harmonize_workflow {
    String connected
    File inputData
    File rnaExpr
    File cnaExpr
    String dataType
    String? analysisDir

  call pgdac_harmonize {
    input:
      inputData=inputData,
      rnaExpr=rnaExpr,
      cnaExpr=cnaExpr,
      analysisDir=analysisDir,
      connected=connected,
      type=dataType
  }
}
