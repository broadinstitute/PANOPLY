task pgdac_rna_protein_correlation_report {
  File tarball
  String label
  String type = "proteome"
  String tmpDir = "tmp"
  Float fdr = 0.05

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /src/rmd-rna-seq-correlation.r ${tarball} ${label} ${type} ${fdr} ${tmpDir}
  }

  output {
    File report = "rna-corr.html"
  }

  runtime {
    docker : "broadcptac/pgdac_cpdb:4"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}


workflow pgdac_rna_protein_correlation_report_workflow {
	call pgdac_rna_protein_correlation_report
}
