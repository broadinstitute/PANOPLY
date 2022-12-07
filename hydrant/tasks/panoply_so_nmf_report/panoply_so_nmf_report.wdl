#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_so_nmf_report {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    File nmf_tar
    File sankey_tar
    String label

    command {
        set -euo pipefail

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/so-nmf-renderRMD.R "${nmf_tar}" "${sankey_tar}" "${label}"
    }

    output {
        File report_out = "${label}_so_nmf_rmd.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_so_nmf_report:latest"
        memory : select_first ([memory, 10]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
        cpu   : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "C. Williams"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_so_nmf_report_workflow {
    call panoply_so_nmf_report
}

