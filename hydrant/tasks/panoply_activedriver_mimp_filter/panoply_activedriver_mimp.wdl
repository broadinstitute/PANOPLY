workflow panoply_activedriver_mimp_filter_workflow {

    call panoply_activedriver_mimp_filter

}

task panoply_activedriver_mimp_filter {
    File mimp_tar
    File activedriver_tar
    String output_prefix
    
    Int? num_cpus
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    command {
        set -euo pipefail

        cp -s ${mimp_tar} .
        cp -s ${activedriver_tar} .

		/usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/panoply_activedriver_mimp_filter.R

		tar -czvf "${output_prefix}_mimp_results_filtered_dir.tar" mimp_results_filtered_dir
    }

    output {
        File tar_out = "${output_prefix}_mimp_results_filtered_dir.tar"
    }

    runtime {
        docker: "broadcptacdev/panoply_activedriver_mimp_filter:latest"
        cpu: "${if defined(num_cpus) then num_cpus else '2'}"
        memory: "${if defined(ram_gb) then ram_gb else '8'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Khoi Pham Munchic"
        email : "proteogenomics@broadinstitute.org"
    }
}
