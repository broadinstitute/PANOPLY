workflow panoply_omicsev_workflow {
    call panoply_omicsev
}


task panoply_omicsev {
	String STANDALONE
    File yaml_file
	String label

	Array[File]? data_files
	File? sample_anno_file
	File? rna_file
    File? panoply_harmonize_tar_file

    String? class_column_name
    String? batch_column_name
	Boolean? data_log_transformed
	Boolean? rna_log_transformed
    Boolean? do_function_prediction
    
    Int? cpu
    Int? memory
    Int? local_disk_gb
    Int? num_preemptions

    command {
	set -euo pipefail

	if [ ${STANDALONE} == "false" ]; then
		mkdir panoply_harmonize_output
		tar -xf ${panoply_harmonize_tar_file} -C panoply_harmonize_output
		output_dir="$(pwd)/panoply_harmonize_output/$(ls panoply_harmonize_output | head -1)/omicsev-data"
		tar_dir="$(pwd)/panoply_harmonize_output"
        tar_name="$(ls panoply_harmonize_output | head -1)"
		home_dir="$(pwd)"

		Rscript /prot/proteomics/Projects/PGDAC/src/omicsev/validate_harmonize_tar.R "$(pwd)/panoply_harmonize_output"

		data_files="$(cat data_file.txt)"
		rna_file="$tar_dir/$tar_name/harmonized-data/rna-matrix.csv"
		sample_anno_file="$tar_dir/$tar_name/harmonized-data/sample-info.csv"
	
    else
    
		output_dir="$(pwd)/omicsev-data"
		tar_dir="$(pwd)"
        tar_name="omicsev-data"
		home_dir="$(pwd)"

		data_files="${sep=',' data_files}"
		rna_file=${default='' rna_file}
		sample_anno_file="${sample_anno_file}"
	fi

	mkdir -p $output_dir
	cd $output_dir


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

	if [ ${STANDALONE} == "false" ]; then
		cp final_output_params.yaml $tar_dir/updated-master-parameter.yaml
	fi


	echo "Preprocessing"
    
    if [ -z $rna_file ]; then

      Rscript \
      /prot/proteomics/Projects/PGDAC/src/omicsev/panoply_omicsev_preprocessing.R \
      --STANDALONE ${STANDALONE} \
      --yaml_file final_output_params.yaml \
      --data_files $data_files \
      --sample_anno_file $sample_anno_file
      
	else
    
      Rscript \
      /prot/proteomics/Projects/PGDAC/src/omicsev/panoply_omicsev_preprocessing.R \
      --STANDALONE ${STANDALONE} \
      --yaml_file final_output_params.yaml \
      --data_files $data_files \
      --sample_anno_file $sample_anno_file \
      --rna_file $rna_file
      
    fi


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

	cd $home_dir

	cp "$output_dir/final_evaluation_report.html" "final_evaluation_report.html"
	mv "final_evaluation_report.html" "omicsev_${label}.html"
    
    cd $tar_dir
    tar -czvf "$home_dir/omicsev_output.tar" $tar_name
    
    }

    runtime {
    	docker: "broadcptacdev/panoply_omicsev:latest"
        memory: "${if defined(memory) then memory else '96'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
        cpu: "${if defined(cpu) then cpu else '6'}"
    }
    output {
        File report = "omicsev_" + label + ".html"
        File outputs = "omicsev_output.tar"
    }
}