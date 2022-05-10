#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_blacksheep {
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    File input_gct
    File master_yaml
    String output_prefix

    String? apply_filtering
    File? identifiers_file
    File? groups_file
    Float? fraction_samples_cutoff
    Float? fdr_value

    command {
        set -euo pipefail

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
        --module blacksheep \
        --master_yaml ${master_yaml} \
        ${"--blacksheep_apply_filtering " + apply_filtering} \
        ${"--blacksheep_identifiers_file " + identifiers_file} \
        ${"--blacksheep_groups_file " + groups_file} \
        ${"--blacksheep_fraction_samples_cutoff " + fraction_samples_cutoff} \
        ${"--blacksheep_fdr_value " + fdr_value}

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/blacksheep_rcode.R "${input_gct}" "final_output_params.yaml"

        if [ "${groups_file}" != "" ]; then
            cp ${groups_file} "blacksheep"
        fi

        tar -czvf "${output_prefix}_blacksheep.tar" blacksheep final_output_params.yaml
    }

    output {
        File tar_out = "${output_prefix}_blacksheep.tar"
    }

    runtime {
        docker : "broadcptac/panoply_blacksheep:1_1"
        memory : select_first ([memory, 10]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
        cpu : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}


workflow panoply_blacksheep_workflow {
    call panoply_blacksheep
        
}