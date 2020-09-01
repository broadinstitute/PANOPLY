task panoply_blacksheep {
    Float? ram_gb
    Int? local_disk_gb
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
        docker : "broadcptacdev/panoply_blacksheep:latest"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karen Christianson"
        email : "karen@broadinstitute.org"
    }
}

workflow panoply_blacksheep_workflow {

    call panoply_blacksheep

}
