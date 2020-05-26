task panoply_immune_analysis {
  File inputData
  String type
  String standalone
  String? analysisDir
  File? groupsFile
  String? subType
  File? params
  Float ?fdr

  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_immune_analysis-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh immune \
                  -i ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -o ${outFile} \
                  ${"-g " + groupsFile} \
                  ${"-m " + subType} \
                  ${"-z " + fdr} \
                  ${"-p " + params};
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh immune \
                  -rna ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -r ${analysisDir} \
                  -o ${outFile} \
                  ${"-g " + groupsFile} \
                  ${"-m " + subType} \
                  ${"-z " + fdr} \
                  ${"-p " + params}
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptacdev/panoply_immune_analysis:70df9e7"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

workflow panoply_immune_analysis_workflow {
    call panoply_immune_analysis 
}
