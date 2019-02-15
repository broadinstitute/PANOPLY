task pgdac_rna_protein_correlation_piped {
  File tarball
  File rnaExpr
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_rna_protein_correlation-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr -i ${tarball} -t ${type} -c ${codeDir} -rna ${rnaExpr} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_rna_protein_correlation:1"
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


task pgdac_rna_protein_correlation_sep {
  String? analysisDir
  File filteredData
  File rnaExpr
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_rna_protein_correlation-output.tar"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"  

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr -f ${filteredData} -t ${type} -c ${codeDir} -d ${dataDir} -rna ${rnaExpr} -r ${analysisDir} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_rna_protein_correlation:1"
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

workflow pgdac_rna_protein_correlation_workflow {
    Boolean isPiped
    File rnaExpr
    String dataType
    File inputTarOrFiltered
    String? analysisDir

    if(isPiped){ 
        call pgdac_rna_protein_correlation_piped{
            input:
                tarball=inputTarOrFiltered,
                type=dataType,
                rnaExpr=rnaExpr
        } 
    }
    if(!isPiped){ 
        call pgdac_rna_protein_correlation_sep{
            input:
                analysisDir=analysisDir,
                filteredData=inputTarOrFiltered,
                rnaExpr=rnaExpr,
                type=dataType
        } 
    }
}
