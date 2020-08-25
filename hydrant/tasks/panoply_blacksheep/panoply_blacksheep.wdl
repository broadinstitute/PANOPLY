task panoply_blacksheep {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    File input_gct
    File master_yaml

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
        ${"--apply_filtering " + apply_filtering} \
        ${"--identifiers_file " + identifiers_file} \
        ${"--groups_file " + groups_file} \
        ${"--fraction_samples_cutoff " + fraction_samples_cutoff} \
        ${"--fdr_value " + fdr_value}

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/blacksheep_rcode.R "${input_gct}" "final_output_params.yaml"

        tar -czvf blacksheep_outlier_analysis.tar.gz blacksheep 
    }

    output {
        File tar_out = "blacksheep_outlier_analysis.tar.gz"
        File final_yaml = "final_output_params.yaml"
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
