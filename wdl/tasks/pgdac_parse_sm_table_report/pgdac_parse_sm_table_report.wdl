task pgdac_parse_sm_table_report {
  File tarball
  String label
  String type
  String tmpDir = "tmp"

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-normalize.r ${tarball} ${label} ${type} ${tmpDir}
  }

  output {
    File report = "norm.html"
  }

  runtime {
    docker : "broadcptac/pgdac_rmd:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}


workflow pgdac_parse_sm_table_report_workflow {
	call pgdac_parse_sm_table_report
}
