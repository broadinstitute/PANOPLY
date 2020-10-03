#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_association_report {
    Int? memory
  	Int? disk_space
  	Int? num_threads

    File input_tar
    File master_yaml
    String label
    String type

    Float? fdr_value

    command {
        set -euo pipefail

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
        --module association_report \
        --master_yaml ${master_yaml} \
        ${"--fdr_value " + fdr_value}

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/rmd_association.r "${input_tar}" "final_output_params.yaml" "${label}" "${type}"
    }

    output {
        File report_out = "${label}_${type}_association_rmd.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_association_report:latest"
        memory: select_first ([memory, 8]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
        preemptible : select_first ([num_threads, 1]) + ""
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_association_report_workflow {

    call panoply_association_report
    
}
