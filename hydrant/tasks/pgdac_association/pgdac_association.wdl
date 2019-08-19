task pgdac_association {
  File inputData
  String type
  String standalone
  String? analysisDir
  File? groupsFile
  String? subType
  File? params

  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_association-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh assoc \
                  -i ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -o ${outFile} \
                  ${"-g " + groupsFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh assoc \
                  -f ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -r ${analysisDir} \
                  -o ${outFile} \
                  -g ${groupsFile} \
                  ${"-m " + subType} \
                  ${"-p " + params}
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_association:1"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email : "rkothadi@broadinstitute.org"
  }
}

workflow pgdac_association_workflow {
  String standalone
  File inputData
  String? analysisDir
  File? groupsFile
  String dataType
    
  call pgdac_association {
    input:
      inputData=inputData,
      type=dataType,
      standalone=standalone,
      analysisDir=analysisDir,
      groupsFile=groupsFile
  }
}
