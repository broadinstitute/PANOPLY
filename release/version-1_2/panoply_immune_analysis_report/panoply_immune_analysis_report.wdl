#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_immune_analysis_report {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    File tar_file
    File yaml_file
    String label

    command {
        set -euo pipefail

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/rmd_immune_analysis.r "${tar_file}" "${yaml_file}" "${label}"
    }

    output {
        File report_out = "${label}_immune_rmd.html"
    }

    runtime {
        docker : "broadcptac/panoply_immune_analysis_report:1_2"
        memory : select_first ([memory, 2]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 10]) + " SSD"
        cpu   : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_immune_analysis_report_workflow {

    call panoply_immune_analysis_report 

}