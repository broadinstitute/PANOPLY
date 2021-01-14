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
    String output_prefix

    String? sm_output

    command {
        set -euo pipefail

        Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
        --module ptm_normalization \
        --master_yaml ${yaml} \
        ${"--spectrum_mill_output " + sm_output}

        Rscript /prot/proteomics/Projects/PGDAC/src/normalize-ptm.R \
          ${proteome_gct} \
          ${ptm_gct} \
          ${default="NULL" output_prefix} \
          ${yaml}
    }

    output {
        File tar_out = "${output_prefix}_blacksheep.tar"
    }

    runtime {
        docker : "broadcptacdev/panoply_blacksheep:latest"
        memory : select_first ([memory, 10]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
        cpu : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "D. R. Mani"
        email : "proteogenomics@broadinstitute.org"
    }
}


workflow panoply_blacksheep_workflow {
    call panoply_blacksheep
        
}