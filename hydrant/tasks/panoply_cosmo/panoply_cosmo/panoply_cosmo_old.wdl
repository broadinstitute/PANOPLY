workflow panoply_cosmo_workflow {
	panoply_cosmo
}

task panoply_cosmo_task {
	String STANDALONE
    File yaml_file

    String? sample_label
	Boolean? run_cosmo
	File? panoply_harmonize_tar
    File? d1_file
    File? r2_file
    File? sample_file

    Int? cpu
    Int? memory
    Int? local_disk_gb
    Int? num_preemptions
    
    Int? cpu
    Int? memory
    Int? local_disk_gb
    Int? num_preemptions




    String d1_file_base = if (STANDALONE == "true") then basename(d1_file, ".gct") else "data"
    String d2_file_base = if (STANDALONE == "true") then basename(d2_file, ".gct") else "rna"
    String sample_file_base = basename(sample_file, ".csv")


	command {
		set -euo pipefail

		echo "Managing parameters"

		Rscript \
		/prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
		--module cosmo \
		--master_yaml ${yaml_file} \
		${if defined(sample_label) then "--cosmo_sample_label " else ""}${sample_label}

		run_cosmo=$(Rscript /prot/proteomics/Projects/PGDAC/src/cosmo/run_cosmo_from_yaml.R final_output_params.yaml)

		if [ $run_cosmo == "FALSE" ]; then
			echo "COSMO not run"
			cp ${
			exit 0







		if [ ${STANDALONE} == "false" ]; then
			mkdir panoply_harmonize_output
			tar -xf ${panoply_harmonize_tar_file} -C panoply_harmonize_output
			output_dir="$(pwd)/panoply_harmonize_output/$(ls panoply_harmonize_output | head -1)/cosmo-data"
			tar_dir="$(pwd)/panoply_harmonize_output"
        	tar_name="$(ls panoply_harmonize_output | head -1)"
			home_dir="$(pwd)"

			Rscript /prot/proteomics/Projects/PGDAC/src/omicsev/validate_harmonize_tar.R "$(pwd)/panoply_harmonize_output"

			d1_file="$(cat data_file.txt)"
			d2_file="$tar_dir/$tar_name/harmonized-data/rna-matrix.csv"
			sample_file="$tar_dir/$tar_name/harmonized-data/sample-info.csv"
	
    	else
    
			output_dir="$(pwd)/cosmo-data"
			tar_dir="$(pwd)"
        	tar_name="cosmo-data"
			home_dir="$(pwd)"

			d1_file=${d1_file}
			d2_file=${d2_file}
			sample_file=${sample_file}
		fi

		mkdir -p $output_dir
		cd $output_dir




        d1_file_preprocessed="./preprocessed_data/d1.tsv"
        d2_file_preprocessed="./preprocessed_data/d2.tsv"
        sample_file_preprocessed="./preprocessed_data/sample-list.tsv"
        sample_label_preprocessed="./preprocessed_data/sample_label.txt"
        
        cosmo_out_dir="./output"
        data_use_dir="./data_use"

        d1_file_use="$data_use_dir/d1.tsv"
        d2_file_use="$data_use_dir/d2.tsv"
        sample_file_use="$data_use_dir/sample-list.tsv"

        method1_out_folder="$cosmo_out_dir/method1_folder"
        method2_out_folder="$cosmo_out_dir/method2_folder"
        final_res_out_folder="$cosmo_out_dir/final_res_folder"
        
		mkdir $cosmo_out_dir
        mkdir $data_use_dir


		










		--d1_file d1_file \
		--d2_file d2_file \
		--sample_file sample_file \
		${if defined(class_column_name) then "--omicsev_class_column_name " else ""}${class_column_name} \
		${if defined(batch_column_name) then "--omicsev_batch_column_name " else ""}${batch_column_name} \
		${if defined(data_log_transformed) then "--omicsev_data_log_transformed " else ""}${data_log_transformed} \
		${if defined(rna_log_transformed) then "--omicsev_rna_log_transformed " else ""}${rna_log_transformed} \
		${if defined(do_function_prediction) then "--omicsev_do_function_prediction " else ""}${do_function_prediction}

		if [ ${STANDALONE} == "false" ]; then
			cp final_output_params.yaml $tar_dir/updated-master-parameter.yaml
		fi





        
        
        Rscript /opt/src/panoply_cosmo_preprocessing.R \
		--STANDALONE ${STANDALONE} \
		--yaml_file 



        ${panoply_harmonize_tar} \
		${yaml_file} \
		${default="none" sample_label}
        
        sample_label=$(cat $sample_label_preprocessed)
		

		R -e \
		"source('/opt/cosmo/tools.R');
		format_input_data('$d1_file_preprocessed', '$d2_file_preprocessed', '$sample_file_preprocessed', out_dir = '$data_use_dir')"

		
		mkdir $method1_out_folder

		R -e \
        "source('/opt/cosmo/method1_function.R');
        d1_file <- '$d1_file_use'; 
        d2_file <- '$d2_file_use'; 
        sample_file <- '$sample_file_use';
        gene_file <- '/opt/cosmo/genes.tsv';
        out_dir <- '$method1_out_folder';
        clinical_attributes <- unlist(strsplit(x='$sample_label',split=','));
        run_2b(d1_file, d2_file, sample_file, gene_file, out_dir=out_dir, clinical_attributes=clinical_attributes)"
		
		mkdir $method2_out_folder

		python /opt/cosmo/method2_function.py \
        -d1 $d1_file_use \
        -d2 $d2_file_use \
        -s $sample_file_use \
        -l $sample_label \
        -o $method2_out_folder

		mkdir $final_res_out_folder

		R -e \
        "source('/opt/cosmo/method1_function.R');
        source('/opt/cosmo/combine_methods.R');
        method1_folder <- '$method1_out_folder';
        method2_folder <- '$method2_out_folder';
        out_dir <- '$final_res_out_folder';
        sample_annotation_file <- '$sample_file_use';
        clinical_attributes <- unlist(strsplit(x='$sample_label',split=','));
        combine_methods(method1_folder, method2_folder, sample_annotation_file, clinical_attributes = clinical_attributes, out_dir = out_dir, prefix = 'cosmo')"

		zip -r "cosmo_output.zip" $out_dir

	}

	runtime {
		docker: "broadcptacdev/panoply_cosmo:latest"
		memory: "${if defined(memory) then memory else '16'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
        cpu: "${if defined(cpu) then cpu else '6'}"
	}

	output {
		File output_zip = "cosmo_output.zip"
	}
    
    parameter_meta {
    	panoply_harmonize_tar: "Output from panoply_harmonize"
        sample_label: "Column header(s) in sample annotation file to use for prediction. This typically includes 'gender'. Column will be excluded if it has any unique levels (e.g. only one female sample, the rest are male). If multiple inputs, separate inputs with a comma (e.g. 'gender,msi'). If no input, columns are pulled from groups in yaml file"
        yaml_file: "Yaml file with PANOPLY set-up."
    }
    
    meta {
    	author: "Stephanie Vartany"
        email: "proteogenomics@broadinstitute.org"
    }
}

task COSMO_task_STANDALONE {

    File? d1_file
    File? d2_file
    File? sample_file
    String? sample_label
    File yaml_file

    Int? cpu
    Int? memory
    Int? local_disk_gb
    Int? num_preemptions
    
    String d1_file_base = basename(d1_file, ".gct")
    String d2_file_base = basename(d2_file, ".gct")
    String sample_file_base = basename(sample_file, ".csv")

	command {
		set -euo pipefail

        
        d1_file_preprocessed="./cosmo_preprocessed_data/${d1_file_base}.tsv"
        d2_file_preprocessed="./cosmo_preprocessed_data/${d2_file_base}.tsv"
        sample_file_preprocessed="./cosmo_preprocessed_data/${sample_file_base}.tsv"
        sample_label_preprocessed="./cosmo_preprocessed_data/sample_label.txt"
        
        out_dir="./output"
        data_use_dir="./data_use"

        d1_file_use="$data_use_dir/${d1_file_base}.tsv"
        d2_file_use="$data_use_dir/${d2_file_base}.tsv"
        sample_file_use="$data_use_dir/${sample_file_base}.tsv"

        method1_out_folder="$out_dir/method1_folder"
        method2_out_folder="$out_dir/method2_folder"
        final_res_out_folder="$out_dir/final_res_folder"
        
		mkdir $out_dir
        mkdir $data_use_dir
        
        
        Rscript /opt/src/panoply_cosmo_preprocessing_standalone.R \
        ${d1_file} \
        ${d2_file} \
        ${sample_file} \
        ${default="none" sample_label} \
		${yaml_file}
        
        sample_label=$(cat $sample_label_preprocessed)
		

		R -e \
		"source('/opt/cosmo/tools.R');
		format_input_data('$d1_file_preprocessed', '$d2_file_preprocessed', '$sample_file_preprocessed', out_dir = '$data_use_dir')"

		
		mkdir $method1_out_folder

		R -e \
        "source('/opt/cosmo/method1_function.R');
        d1_file <- '$d1_file_use'; 
        d2_file <- '$d2_file_use'; 
        sample_file <- '$sample_file_use';
        gene_file <- '/opt/cosmo/genes.tsv';
        out_dir <- '$method1_out_folder';
        clinical_attributes <- unlist(strsplit(x='$sample_label',split=','));
        run_2b(d1_file, d2_file, sample_file, gene_file, out_dir=out_dir, clinical_attributes=clinical_attributes)"
		
		mkdir $method2_out_folder

		python /opt/cosmo/method2_function.py \
        -d1 $d1_file_use \
        -d2 $d2_file_use \
        -s $sample_file_use \
        -l $sample_label \
        -o $method2_out_folder

		mkdir $final_res_out_folder

		R -e \
        "source('/opt/cosmo/method1_function.R');
        source('/opt/cosmo/combine_methods.R');
        method1_folder <- '$method1_out_folder';
        method2_folder <- '$method2_out_folder';
        out_dir <- '$final_res_out_folder';
        sample_annotation_file <- '$sample_file_use';
        clinical_attributes <- unlist(strsplit(x='$sample_label',split=','));
        combine_methods(method1_folder, method2_folder, sample_annotation_file, clinical_attributes = clinical_attributes, out_dir = out_dir, prefix = 'cosmo')"

		zip -r "cosmo_output.zip" $out_dir

	}

	runtime {
		docker: "broadcptacdev/panoply_cosmo:latest"
		memory: "${if defined(memory) then memory else '16'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
        cpu: "${if defined(cpu) then cpu else '6'}"
	}

	output {
		File output_zip = "cosmo_output.zip"
	}
    
    parameter_meta {
    	d1_file: "First file to analyze for mislabeling (typically at the protein-level). Must be a GCT."
        d2_file: "Second file to analyze for mislabeling (typically transcriptome). Must be a GCT."
        sample_file: "Sample annotation file. Must be .csv file."
        sample_label: "Column header(s) in sample annotation file to use for prediction. This typically includes 'gender'. Column will be excluded if it has any unique levels (e.g. only one female sample, the rest are male). If multiple inputs, separate inputs with a comma (e.g. 'gender,msi'). If no input, columns are pulled from groups in yaml file."
        yaml_file: "Yaml file with PANOPLY set-up."
    }
    
    meta {
    	author: "Stephanie Vartany"
        email: "proteogenomics@broadinstitute.org"
    }
}