version 1.0

workflow panoply_omicsev_workflow {
    call panoply_omicsev

}


task panoply_omicsev {
	input {
		String STANDALONE
	  File yaml_file
		String label

		Array[File]? data_files
		File? sample_anno_file
		File? rna_file
	  File? panoply_harmonize_tar_file
	  String? ome_type

	  String? class_column_name
	  String? batch_column_name
		Boolean? data_log_transformed
		Boolean? rna_log_transformed
	  Boolean? do_function_prediction

	  Int? memory
	  Int? disk_space
	  Int? num_threads
	  Int? num_preemptions
	}

  command {
  	set -euo pipefail
  
  	if [ ${STANDALONE} == "false" ]; then
  	
  	  tar -xf ${select_first([panoply_harmonize_tar_file, ""])}
  	  
  	  # get the root directory of the tar file
  	  tar -tf ${select_first([panoply_harmonize_tar_file, ""])} > all_files_in_tar.txt
  	  tar_dir=$(pwd)/$(basename $(head -n 1 all_files_in_tar.txt))
  	  
  	  Rscript /prot/proteomics/Projects/PGDAC/src/omicsev/validate_harmonize_tar.R $tar_dir ${select_first([ome_type, ""])}
  	  
  	  data_files="$tar_dir/harmonized-data/${select_first([ome_type, ""])}-matrix.csv"
  		rna_file="$tar_dir/harmonized-data/rna-matrix.csv"
  		sample_anno_file="$tar_dir/harmonized-data/sample-info.csv"
  	
    else
      data_files="${if defined(data_files) then sep(',', select_first([data_files])) else ""}"
  		rna_file=${select_first([rna_file, ''])}
  		sample_anno_file="${select_first([sample_anno_file, ""])}"
  	fi
  
    output_dir="$(pwd)/omicsev-data"
  	mkdir -p $output_dir
  	cd $output_dir
  
  
  	echo "Managing parameters"
  
  	Rscript \
    	/prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    	--module omicsev \
    	--master_yaml ${yaml_file} \
    	${"--omicsev_class_column_name " + class_column_name} \
    	${"--omicsev_batch_column_name " + batch_column_name} \
    	${if defined(data_log_transformed) then "--omicsev_data_log_transformed " + select_first([data_log_transformed]) else ""} \
    	${if defined(rna_log_transformed) then "--omicsev_rna_log_transformed " + select_first([rna_log_transformed]) else ""} \
    	${if defined(do_function_prediction) then "--omicsev_do_function_prediction " + select_first([do_function_prediction]) else ""}
  
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
        ${select_first([num_threads, 6])} \
        protein \
        x2.tsv \
        $(cat do_function_prediction.txt) \
        "./"
  
  	cd ..
  
  	cp "$output_dir/final_evaluation_report.html" "final_evaluation_report.html"
  	mv "final_evaluation_report.html" "omicsev_${label}.html"
      
    tar -czvf "omicsev_output.tar" $(basename $output_dir)
      
    }

    runtime {
    	docker: "broadcptacdev/panoply_omicsev:latest"
        memory: "${select_first([memory, 96])}GB"
        disks : "local-disk ${select_first([disk_space, 10])} HDD"
        preemptible : select_first([num_preemptions, 0])
        cpu : select_first([num_threads, 6])
    }
    output {
        File report = "omicsev_" + label + ".html"
        File outputs = "omicsev_output.tar"
    }
}