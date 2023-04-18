task panoply_activedriver {
    Int? num_cpus
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    File phospho_file
    File mutation_file
    File fasta_file
    File master_yaml
    String output_prefix

	String? search_engine
	String? group_by_column
    String? residues
	String? mutation_AA_change_colname
	String? sample_id_col
	String? transcript_id_col
    String? other_transcript_col
    String? region_pval_adj_method

    command {
        set -euo pipefail

		/usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/panoply_activedriver.R "${phospho_file}" "${mutation_file}" "${fasta_file}" "${master_yaml}" "${num_cpus}"

		tar -czvf "${output_prefix}_activedriver_output.tar" activedriver_results_dir
    }

    output {
        File tar_out = "${output_prefix}_activedriver_output.tar"
    }

    runtime {
        docker: "broadcptacdev/panoply_activedriver:latest"
        cpu: "${if defined(num_cpus) then num_cpus else '8'}"
        memory: "${if defined(ram_gb) then ram_gb else '32'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Khoi Pham Munchic"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_activedriver_workflow {

    call panoply_activedriver 

}
