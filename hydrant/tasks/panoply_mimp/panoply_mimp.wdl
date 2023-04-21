task panoply_mimp {
    Int? num_cpus
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    File mutation_file
    File phospho_file
    File fasta_file
    File ids_file
    File master_yaml
    String output_prefix
    File? kinase_model

    File? groups_file_path
	String? search_engine
	String? phosphosite_col
    String? protein_id_col
    String? protein_id_type
	String? mutation_AA_change_colname
	String? mutation_type_col
	String? sample_id_col 
	String? transcript_id_col

    command {
        set -euo pipefail

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
        --module mimp \
        --master_yaml ${master_yaml} \
        ${"--mimp_groups_file_path " + groups_file_path} \
        ${"--mimp_search_engine " + search_engine} \
        ${"--mimp_phosphosite_col " + phosphosite_col} \
        ${"--mimp_protein_id_col " + protein_id_col} \
        ${"--mimp_protein_id_type " + protein_id_type} \
        ${"--mimp_mutation_AA_change_colname " + mutation_AA_change_colname} \
        ${"--mimp_mutation_type_col " + mutation_type_col} \
        ${"--mimp_sample_id_col " + sample_id_col} \
        ${"--mimp_transcript_id_col " + transcript_id_col}

        /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/panoply_mimp.R "${mutation_file}" "${phospho_file}" "${fasta_file}" "${ids_file}" "final_output_params.yaml" "${kinase_model}"

		tar -czvf "${output_prefix}_mimp_output.tar" mimp_results_dir

    }

    output {
        File tar_out = "${output_prefix}_mimp_output.tar"
    }

    runtime {
        docker : "broadcptacdev/activedriver_mimp:latest"
        cpu: "${if defined(num_cpus) then num_cpus else '8'}"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karen Christianson"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_mimp_workflow {

    call panoply_mimp 

}
