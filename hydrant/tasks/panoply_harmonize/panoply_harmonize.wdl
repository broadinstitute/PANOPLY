task panoply_harmonize {
  File inputData
  File rnaExpr
  File cnaExpr
  String type
  String standalone
  String? analysisDir
  String? subType
  File? params

  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"
  String outFile = "panoply_harmonize-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    if [[ ${standalone} = false ]]; then
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
    docker : "broadcptacdev/panoply_harmonize:latest"
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

workflow panoply_harmonize_workflow {
    String standalone
    File inputData
    File rnaExpr
    File cnaExpr
    String dataType
    String? analysisDir

  call panoply_harmonize {
    input:
      inputData=inputData,
      rnaExpr=rnaExpr,
      cnaExpr=cnaExpr,
      analysisDir=analysisDir,
      standalone=standalone,
      type=dataType
  }
}
