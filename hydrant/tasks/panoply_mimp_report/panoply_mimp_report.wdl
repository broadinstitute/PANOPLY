version 1.0

task panoply_mimp_report {
    input {
      Float? memory
      Int? disk_space
      Int? num_preemptions

      File tar_file
      String output_prefix
    }

    command {
        set -euo pipefail

        Rscript /prot/proteomics/Projects/PGDAC/src/rmd_mimp.R "${tar_file}" "${output_prefix}"
    }

    output {
        File report_out = "${output_prefix}_mimp_rmd.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_mimp_report:latest"
        memory: "${if defined(memory) then memory else '2'}GB"
        disks : "local-disk ${if defined(disk_space) then disk_space else '10'} HDD"
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
