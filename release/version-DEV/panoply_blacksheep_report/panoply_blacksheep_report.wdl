#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_blacksheep_report {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    File input_tar
    String output_prefix
    String type

    command {
        set -euo pipefail

        /usr/bin/Rscript /home/pgdac/src/rmd_blacksheep.R "${input_tar}" "${output_prefix}" "${type}"
    }

    output {
        File report_out = "${output_prefix}_blacksheep_rmd.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_blacksheep_report:DEV"
        memory : select_first ([memory, 10]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
        cpu : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_blacksheep_report_workflow {

    call panoply_blacksheep_report

}

