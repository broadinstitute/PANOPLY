task pgdac_cons_clust {
  File tarball   # output from pgdac_harmonize or pgdac_normalize_ms_data
  String type
  File? groupsFile
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_cluster-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    echo ${type}
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh cluster -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-m " + subType} ${"-p " + params} ${"-g " + groupsFile}
   
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_cons_clust:2"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani, Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}


workflow pgdac_cons_clust_workflow {
	call pgdac_cons_clust
}
