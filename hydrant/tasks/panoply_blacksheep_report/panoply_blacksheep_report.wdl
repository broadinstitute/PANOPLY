task panoply_blacksheep_report {
    Int? memory
  	Int? disk_space
  	Int? num_threads

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
        memory : select_first ([memory, 10]) + "GB"
    	disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    	cpu : select_first ([num_threads, 1]) + ""
    }

    meta {
        author : "Karen Christianson"
        email : "karenchristianson@broadinstitute.org"
    }
}

workflow panoply_blacksheep_report_workflow {

    call panoply_blacksheep_report

}
