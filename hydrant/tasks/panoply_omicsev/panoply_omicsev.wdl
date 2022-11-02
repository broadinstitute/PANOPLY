workflow panoply_omicsev_workflow {
    call OmicsEV_task {}
}


task OmicsEV_task {
	String STANDALONE
    File yaml_file

	Array[File]? data_files
	File? sample_anno_file
	File? rna_file
    File? panoply_harmonize_tar_file

    String? class_column_name
    String? batch_column_name
	Boolean? data_log_transformed
	Boolean? rna_log_transformed
    Boolean? do_function_prediction

    String output_dir = "omicsev_output"
    
    Int? cpu
    Int? memory
    Int? local_disk_gb
    Int? num_preemptions

    command {
	set -euo pipefail
    
    mkdir ${output_dir}
    
    cd ${output_dir}

	echo "Managing parameters"

	Rscript \
	/prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
	--module omicsev \
	--master_yaml ${yaml_file} \
	${if defined(class_column_name) then "--omicsev_class_column_name " else ""}${class_column_name} \
	${if defined(batch_column_name) then "--omicsev_batch_column_name " else ""}${batch_column_name} \
	${if defined(data_log_transformed) then "--omicsev_data_log_transformed " else ""}${data_log_transformed} \
	${if defined(rna_log_transformed) then "--omicsev_rna_log_transformed " else ""}${rna_log_transformed} \
	${if defined(do_function_prediction) then "--omicsev_do_function_prediction " else ""}${do_function_prediction}

	echo "Preprocessing"
    
    Rscript \
    /prot/proteomics/Projects/PGDAC/src/omicsev/panoply_omicsev_preprocessing.R \
	--STANDALONE ${STANDALONE} \
	--yaml_file final_output_params.yaml \
	${if defined(data_files) then "--data_files " else ""}${sep=',' data_files} \
	${if defined(sample_anno_file) then "--sample_anno_file " else ""}${sample_anno_file} \
	${if defined(rna_file) then "--rna_file " else ""}${rna_file} \
	${if defined(panoply_harmonize_tar_file) then "--harmonize_tar_file " else ""}${panoply_harmonize_tar_file}

	echo "Running OmicsEV"

    Rscript \
	/prot/proteomics/Projects/PGDAC/src/omicsev/panoply_run_OmicsEV.R \
    dataset \
    sample_list.tsv \
    ${default=6 cpu} \
    protein \
	x2.tsv \
    $(cat do_function_prediction.txt) \
    "./"
    
    cd ..
    
    zip -r "${output_dir}.zip" ${output_dir}
    
    }

    runtime {
    	docker: "broadcptacdev/panoply_omicsev:latest"
        memory: "${if defined(memory) then memory else '96'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
        cpu: "${if defined(cpu) then cpu else '6'}"
    }
    output {
        File html_report = "${output_dir}/final_evaluation_report.html"
        File omicsev_output_zip = "${output_dir}.zip"
    }
}