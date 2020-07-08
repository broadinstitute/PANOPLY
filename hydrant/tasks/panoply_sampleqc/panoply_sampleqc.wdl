task panoply_sampleqc {
  File tarball   # output from panoply_harmonize
  String type
  String? subType
  File yaml
  Float? corThreshold
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_sampleqc-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module sample_qc \
    --master_yaml ${yaml} \
    ${"--cor_threshold " + corThreshold} 
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh sampleQC -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-m " + subType} -p "config-custom.r"
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/panoply_sampleqc:dev"
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


workflow panoply_sampleqc_workflow {
	call panoply_sampleqc
}
