task panoply_rna_protein_correlation {
  File inputData
  File rnaExpr
  String type
  String standalone
  String? subType
  File? params
  String? analysisDir
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_rna_protein_correlation-output.tar"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr \
                  -i ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -rna ${rnaExpr} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr \
                  -f ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -d ${dataDir} \
                  -rna ${rnaExpr} \
                  -r ${analysisDir} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/panoply_rna_protein_correlation:1"
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

workflow panoply_rna_protein_correlation_workflow {
  File rnaExpr
  String dataType
  File inputData
  String standalone
  String? analysisDir

  call panoply_rna_protein_correlation {
    input:
      inputData=inputData,
      type=dataType,
      rnaExpr=rnaExpr,
      analysisDir=analysisDir,
      standalone=standalone
  } 
}
