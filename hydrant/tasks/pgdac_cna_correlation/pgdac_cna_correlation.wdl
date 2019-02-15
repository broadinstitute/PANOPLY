task pgdac_cna_correlation {
  File tarball   # output from pgdac_cna_setup
  String type
  String? subType
  File? params
  String outFile = "pgdac_cna_correlation-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CNAcorr -i ${tarball} -t ${type} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_cna_setup:1"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}


workflow pgdac_cna_correlation_workflow {
	call pgdac_cna_correlation
}
