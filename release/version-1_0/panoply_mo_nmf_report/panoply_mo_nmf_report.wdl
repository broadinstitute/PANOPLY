#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_mo_nmf_report {

    File tarball
    String label

    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    command {
        set -euo pipefail
	Rscript /home/pgdac/src/rmd-mo-nmf.R ${tarball} ${label}
    }

    output {
	File report = "mo-nmf-report-" + label + ".html"
    }

    runtime {
        docker : "broadcptac/panoply_mo_nmf_report:1_0"
        memory: select_first ([memory, 2]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 10]) + " SSD"
        cpu : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Karsten Krug"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_mo_nmf_report_workflow {

    call panoply_mo_nmf_report
}
