#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_mo_nmf_report {

    File tarball
    String label

    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    command {
        set -euo pipefail
	Rscript /home/pgdac/src/rmd-mo-nmf.R ${tarball} ${label}
    }

    output {
	File report = "mo-nmf-report-" + label + ".html"
    }

    runtime {
        docker : "broadcptacdev/panoply_mo_nmf_report:latest"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karsten Krug"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_mo_nmf_report_workflow {

    call panoply_mo_nmf_report
}
