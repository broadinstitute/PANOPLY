#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cons_clust_report {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    File tar_file
    File yaml_file
    String label
    String type

    command {
        set -euo pipefail

        Rscript /prot/proteomics/Projects/PGDAC/src/rmd-cons-clust.r "${tar_file}" "${yaml_file}" "${label}" "${type}"
    }

    output {
        File report_out = "${label}_${type}_cons_clust_rmd.html"
    }

    runtime {
        docker : "broadcptac/panoply_cons_clust_report:1_5"
        memory : select_first ([memory, 10]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
        cpu   : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_cons_clust_report_workflow {

    call panoply_cons_clust_report
}

