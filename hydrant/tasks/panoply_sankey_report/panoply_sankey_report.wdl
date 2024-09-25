#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_sankey_report {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    String annot_of_comparison

    String label
    File sankey_tar

    command {
        set -euo pipefail

        Rscript /prot/proteomics/Projects/PGDAC/src/sankey-renderRMD.R "${annot_of_comparison}" "${sankey_tar}" "${label}"
    }

    output {
        File report_out = "${label}_sankey_rmd.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_sankey_report:latest"
        memory : select_first ([memory, 10]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
        cpu   : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "C.M. Williams"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_sankey_report_workflow {
    call panoply_sankey_report
}