task pgdac_cmap_input {
  File tarball   # output from pgdac_cna_correlation
  String? cmap_group
  String? cmap_type
  String? cmap_log
  Int? cmap_permutations
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_cmapsetup-output.tar"
  String outGmtFile = "cmap-trans-genesets.gmt"
  
  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions
  
  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CMAPsetup -i ${tarball} -c ${codeDir} -o ${outFile} ${"-CMAPgroup " + cmap_group} ${"-CMAPtype " + cmap_type} ${"-CMAPnperm " + cmap_permutations} ${"-CMAPlog " + cmap_log}
  }
  
  output {
    File outputs = "${outFile}"
    File genesets = "${outGmtFile}"
    Array[File] permuted_genesets = glob ("*-permuted-genes-*.gmt")
  }

  runtime {
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

workflow pgdac_cmap_input_workflow {
	call pgdac_cmap_input
}
