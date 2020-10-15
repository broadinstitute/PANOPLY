#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cons_clust_report {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    File tar_file
    File yaml_file
    String label
    String type

    command {
        set -euo pipefail

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/rmd_cons_clust.r "${tar_file}" "${yaml_file}" "${label}" "${type}"
    }

    output {
        File report_out = "${label}_${type}_cons_clust_rmd.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_cons_clust_report:latest"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_cons_clust_report_workflow {

    call panoply_cons_clust_report
}
