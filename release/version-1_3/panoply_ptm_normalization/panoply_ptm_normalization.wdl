#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_ptm_normalization {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    File proteome_gct
    File ptm_gct
    File yaml
    
    String? output_prefix = basename (ptm_gct, ".gct")

    String? accession_number_col
    String? accession_numbers_col
    String? accession_numbers_sep
    String? score_col

    command {
        set -euo pipefail

        Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
        --module ptm_normalization \
        --master_yaml ${yaml} \
        ${"--accession_number_colname " + accession_number_col} \
        ${"--accession_numbers_colname " + accession_numbers_col} \
        ${"--accession_numbers_separator " + accession_numbers_sep} \
        ${"--score_colname" + score_col}


        Rscript /prot/proteomics/Projects/PGDAC/src/normalize-ptm.R \
          ${proteome_gct} \
          ${ptm_gct} \
          ${output_prefix} \
          "final_output_params.yaml"
    }

    output {
        File tar_out = "${output_prefix}-proteome-relative-norm.gct"
    }

    runtime {
        docker : "broadcptac/panoply_ptm_normalization:1_3"
        memory : select_first ([memory, 48]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 60]) + " SSD"
        cpu : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "D. R. Mani"
        email : "proteogenomics@broadinstitute.org"
    }
}


workflow panoply_ptm_normalization_workflow {
    call panoply_ptm_normalization
}
