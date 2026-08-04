
task panoply_mimp_report {
      Float? memory
      Int? disk_space
      Int? num_preemptions

      File tar_file
      String output_prefix

    command {
        set -euo pipefail

        Rscript /prot/proteomics/Projects/PGDAC/src/rmd_mimp.R "${tar_file}" "${output_prefix}"
    }

    output {
        File report_out = "${output_prefix}_mimp_rmd.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_mimp_report:latest"
        memory: "${select_first([memory, 2])}GB"
        disks : "local-disk ${select_first([disk_space, 10])} HDD"
        preemptible : select_first([num_preemptions, 0])
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_mimp_report_workflow {

    call panoply_mimp_report

}
