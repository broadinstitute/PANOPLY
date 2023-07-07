#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_sankey {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    Array[File]+ annot_files          # annotation file(s)
    Array[String]+ annot_file_labels  # datatype(s) / label(s) for the provided annotation file(s)

    String annot_column               # annotation column which should be used for sankey comparisons (e.g. "NMF.consensus")
    String annot_prefix               # prefix to prepent to annotation values (e.g. "C" -> C1 C2 C2, instead of 1 2 3)

    File? annot_file_primary          # annotation file which should be centered / highlighted in comparisons
    String? annot_label_primary       # corresponding label

    String label

    command {
        set -euo pipefail

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/sankey.r "${label}" "${sep="," annot_files}" "${sep="," annot_file_labels}" "${annot_column}" "${annot_prefix}" "${annot_file_primary}" "${annot_label_primary}"

        tar -czvf "${label}_sankey_diagrams.tar" sankey-*html
    }

    output {
        File tar_out = "${label}_sankey_diagrams.tar"
    }

    runtime {
        docker : "broadcptacdev/panoply_sankey:test"
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

workflow panoply_sankey_workflow {
    call panoply_sankey
}