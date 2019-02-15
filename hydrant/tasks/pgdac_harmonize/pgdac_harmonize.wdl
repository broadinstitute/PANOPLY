task pgdac_harmonize_piped {
  File tarball
  File rnaExpr
  File cnaExpr
  String type
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
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize -i ${tarball} -t ${type} -c ${codeDir} -d ${dataDir} -rna ${rnaExpr} -cna ${cnaExpr} -o ${outFile} ${"-m " + subType} ${"-p " + params}
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

task pgdac_harmonize_sep {
  String? analysisDir
  File filteredData
  File rnaExpr
  File cnaExpr
  String type
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
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize -f ${filteredData} -r ${analysisDir} -t ${type} -c ${codeDir} -d ${dataDir} -rna ${rnaExpr} -cna ${cnaExpr} -o ${outFile} ${"-m " + subType} ${"-p " + params}
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
    Boolean isPiped
    File inputTarOrFiltered
    File rnaExpr
    File cnaExpr
    String dataType
    String? analysisDir

    if(isPiped){
        call pgdac_harmonize_piped{
            input:
                tarball=inputTarOrFiltered,
                rnaExpr=rnaExpr,
                cnaExpr=cnaExpr,
                type=dataType
        }
    }
    if(!isPiped){
        call pgdac_harmonize_sep{
            input:
                analysisDir=analysisDir,
                filteredData=inputTarOrFiltered,
                rnaExpr=rnaExpr,
                cnaExpr=cnaExpr,
                type=dataType
        }
    }
}
