task panoply_mimp_report {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    File tar_file
    String output_prefix

    command {
        set -euo pipefail

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/rmd_mimp.R "${tar_file}" "${output_prefix}"
    }

    output {
        File report_out = "${output_prefix}_mimp_rmd.html"
    }

    runtime {
        docker : "broadcptac/panoply_mimp_report:1_4"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_mimp_report_workflow {

    call panoply_mimp_report
}
