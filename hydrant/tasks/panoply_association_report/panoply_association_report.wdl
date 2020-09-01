task panoply_association_report {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

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
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karen Christianson"
        email : "karenchristianson@broadinstitute.org"
    }
}

workflow panoply_association_report_workflow {

    call panoply_association_report
    
}
