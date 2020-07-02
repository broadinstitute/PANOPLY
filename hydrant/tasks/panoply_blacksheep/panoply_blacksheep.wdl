task panoply_blacksheep {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    File input_gct
    File master_yaml

    File? identifiers_file
    File? groups_file
    Float? fraction_cutoff
    Float? fdr_cutoff

    command {
        set -euo pipefail

	# Myranda's parameter integration script: input master_yaml, output yaml_file

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/blacksheep_rcode.R "${input_gct}" "${yaml_file}"

        tar -czvf blacksheep_outlier_analysis.tar.gz blacksheep 
    }

    output {
        File tar_out = "blacksheep_outlier_analysis.tar.gz"
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
