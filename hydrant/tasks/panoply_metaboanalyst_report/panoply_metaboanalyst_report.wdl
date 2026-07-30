#
# Copyright (c) 2024 The Broad Institute, Inc. All rights reserved.
#
version 1.1

task panoply_metaboanalyst_report {
    input {
      String label

      # inputs from panoply_metaboanalyst
      File metaboanalyst_results                ## Rdata file containing results of metaboanalyst()

      Int? memory
      Int? disk_space
      Int? num_threads
      Int? num_preemptions
    }

    command {
        set -euo pipefail
        Rscript /prot/proteomics/Projects/PGDAC/src/metaboanalyst-renderRMD.R ${"-n " + metaboanalyst_results} ${"-x " + label}
    }

    output {
        File report = label + "_metaboanalyst_report.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_metaboanalyst_report:DEV"
        memory: select_first ([memory, 16]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 10]) + " SSD"
        cpu : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "C.M. Williams"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_metaboanalyst_report_workflow {
    call panoply_metaboanalyst_report

}
