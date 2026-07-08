#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
version 1.0

task panoply_clumps_ptm_report {
    input {
      String label

      # output(s) from panoply_clumps_ptm_postprocess
      Array[File]+ postprocess_results                ## Tar file(s) with all results from ClumpsPTM

      File? mapping_params                            ## (optional) mapping parameter file, if it's available

      Int? memory
      Int? disk_space
      Int? num_threads
      Int? num_preemptions
    }

    command {
        set -euo pipefail
        Rscript /prot/proteomics/Projects/PGDAC/src/clumps_ptm-renderRMD.R -i ${sep="," postprocess_results} ${"-m " + mapping_params} ${"-x " + label}
    }

    output {
        File report = label + "_clumps_ptm_report.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_clumps_ptm_report:latest"
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

workflow panoply_clumps_ptm_report_workflow {
    call panoply_clumps_ptm_report

}
