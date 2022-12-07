#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_so_nmf_sankey {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    File tar_file
    String label

    command {
        set -euo pipefail

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/so-nmf-sankey.r "${tar_file}"

        tar -czvf "${label}_sankey_diagrams.tar" sankey-*html
    }

    output {
        File tar_out = "${label}_sankey_diagrams.tar"
    }

    runtime {
        docker : "broadcptacdev/panoply_so_nmf_sankey:latest"
        memory : select_first ([memory, 10]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
        cpu   : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "C. Williams"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_so_nmf_sankey_workflow {

    call panoply_so_nmf_sankey
}

