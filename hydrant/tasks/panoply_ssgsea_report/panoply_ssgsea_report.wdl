task panoply_ssgsea_report {

  File tarball
  File cfg_yaml
  String label

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-ssgsea.R -t ${tarball} -l ${label} -y ${cfg_yaml} -f NA -n NA -p NA -z /home/pgdac/src/
  }

  output {
    File report = "report_" + label + ".html"
  }

  runtime {
    docker : "broadcptacdev/panoply_ssgsea_report:latest"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}

workflow panoply_ssgsea_report_workflow {
	call panoply_ssgsea_report
}
