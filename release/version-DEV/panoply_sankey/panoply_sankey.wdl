#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_sankey {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    Array[File]+ annot_files            # annotation file(s)
    Array[String]+ annot_file_labels    # datatype(s) / label(s) for the provided annotation file(s)
    
    File? annot_file_primary            # annotation file which should be centered / highlighted in comparisons (optional)
    String? annot_label_primary         # corresponding label (optional)

    String? id_column                   # id column for identifying entries (e.g. "Sample.ID"); uses rownames if not provided
    String annot_column                 # annotation column for sankey comparisons (e.g. "NMF.consensus")
    String? annot_prefix                # prefix to prepent to annotation values (e.g. "C" -> C1 C2 C3, instead of 1 2 3)
    String? color_str                   # string of colors (comma separated) to use for sankey diagrams (e.g. '#AA0000,#0000AA' will create a gradient from red to blue)

    String label

    command {
        set -euo pipefail

        Rscript /prot/proteomics/Projects/PGDAC/src/sankey.r -x ${label} -f ${sep="," annot_files} -l ${sep="," annot_file_labels} ${"-j " + annot_file_primary} ${"-m " + annot_label_primary} ${"-i " + id_column} -a ${annot_column} ${"-p " + annot_prefix} ${"-c '" + color_str + "'"}
        tar -czvf "${label}_sankey_diagrams.tar" sankey-*html sankey-*pdf sankey_labels.txt # tar all HTMLs and PDFs, and the TXT with all file-labels
    }

    output {
        File tar_out = "${label}_sankey_diagrams.tar"
    }

    runtime {
        docker : "broadcptacdev/panoply_sankey:DEV"
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