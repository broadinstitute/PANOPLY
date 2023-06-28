#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
task panoply_nmf_balance_omes {

    String label
    Array[File]+ ome_gcts
    Array[String]+ ome_labels

    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    Float? tol
    Float? var
    String? zscore_mode

    command {
        set -euo pipefail

        ## run balanace filter
        Rscript /home/pgdac/src/filter-gcts-to-balance-omes.R -f ${sep="," ome_gcts} -l ${sep="," ome_labels} -t ${default="0.01" tol} -v ${default="0.9" var} -z ${default="rowcol" zscore_mode}
       
    }

    output {
       Array[File]+ omes_balanced=ome_gcts
       File pdf="balance-omes.pdf"
    }

    runtime {
        docker : "broadcptacdev/panoply_nmf_balance_omes:test"
        memory: select_first ([memory, 16]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 10]) + " SSD"
        cpu   : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Karsten Krug"
        email : "proteogenomics@broadinstitute.org"
    }
}

## workflow
workflow panoply_nmf_balance_omes_workflow {
    call panoply_nmf_balance_omes
}

