task panoply_blacksheep_report {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    File input_tar
    String output_prefix

    command {
        set -euo pipefail

        /usr/bin/Rscript /home/pgdac/src/rmd_blacksheep.R "${input_tar}" "${output_prefix}"
    }

    output {
        File report_out = "${output_prefix}_blacksheep_rmd.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_blacksheep_report:latest"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karen Christianson"
        email : "karenchristianson@broadinstitute.org"
    }
}

workflow panoply_blacksheep_report_workflow {

    call panoply_blacksheep_report

}
