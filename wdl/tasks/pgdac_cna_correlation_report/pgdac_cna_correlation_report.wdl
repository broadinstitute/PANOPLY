task pgdac_cna_correlation_report {
  File tarball
  String label
  String type
  String tmpDir
  Float fdr

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-cna-analysis.r ${tarball} ${label} ${type} ${fdr} ${tmpDir}
  }

  output {
    File report = "cna-analysis_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/pgdac_rmd:3"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}

workflow pgdac_cna_correlation_report_workflow {
	call pgdac_cna_correlation_report
}
