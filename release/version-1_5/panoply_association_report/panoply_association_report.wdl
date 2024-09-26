#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_association_report {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    Array[File] ssgsea_assoc_tars
    File master_yaml
    String label
    String type

    Float? fdr_value

    command {
        set -euo pipefail

        Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
        --module association_report \
        --master_yaml ${master_yaml} \
        ${"--fdr_value " + fdr_value}

        # compile tars into useable format, equivalent to scatter_processing() in panoply_downloads
        bash /prot/proteomics/Projects/PGDAC/src/compile_tars.sh ${sep=' ' ssgsea_assoc_tars}

        Rscript /prot/proteomics/Projects/PGDAC/src/rmd_association.r "ssgsea_assoc.tar" "final_output_params.yaml" "${label}" "${type}"
    }

    output {
        File report_out = "${label}_${type}_association_rmd.html"
    }

    runtime {
        docker : "broadcptac/panoply_association_report:1_5"
        memory: select_first ([memory, 8]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
        cpu : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_association_report_workflow {

    call panoply_association_report
    
}
