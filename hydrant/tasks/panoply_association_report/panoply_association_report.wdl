task panoply_association_report {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    File input_tar
    File master_yaml

    command {
        set -euo pipefail

        /usr/bin/Rscript /home/pgdac/src/rmd_association.r "${input_tar}" "${master_yaml}"
    }

    output {
        File report_out = "rmd_association.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_association_report:latest"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karen Christianson"
        email : "karenchristianson@broadinstitute.org"
    }
}

workflow panoply_association_report_workflow {

    call panoply_association_report
    
}
