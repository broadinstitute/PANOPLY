workflow panoply_cosmo_workflow {
	String STANDALONE
  File yaml_file
	Boolean run_cosmo = true

	call panoply_cosmo_manage_parameters {
		input:
			yaml_file = yaml_file,
			run_cosmo = run_cosmo
	}

	call panoply_cosmo {
		input:
			yaml_file = panoply_cosmo_manage_parameters.yaml_file_final,
			run_cosmo = panoply_cosmo_manage_parameters.run_cosmo_final,
			STANDALONE = STANDALONE
	}
}

task panoply_cosmo_manage_parameters {
	File yaml_file
	Boolean run_cosmo
	String? sample_label

	Int? cpu
  Int? memory
  Int? local_disk_gb
  Int? num_preemptions

	command {
		set -euo pipefail

		echo "Managing parameters"

		Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
			--module cosmo \
			--master_yaml ${yaml_file} \
			--cosmo_run_cosmo ${run_cosmo} \
			${if defined(sample_label) then "--cosmo_sample_label " + "'${sample_label}'" else ""}


		R -e "cat(yaml::read_yaml('final_output_params.yaml')[['cosmo.params']][['run_cosmo']], file = 'run_cosmo.txt', sep = '\n')"

	}

	runtime {
		docker: "broadcptacdev/panoply_cosmo:latest"
		memory: "${if defined(memory) then memory else '16'}GB"
    disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu: "${if defined(cpu) then cpu else '6'}"
	}

	output {
		Boolean run_cosmo_final = read_boolean("run_cosmo.txt")
		File yaml_file_final = "final_output_params.yaml"
	}
}


task panoply_cosmo {
	String STANDALONE
	File yaml_file
	Boolean run_cosmo

	File? panoply_harmonize_tar
  File? d1_file
  File? d2_file
  File? sample_file

	Int? cpu
  Int? memory
  Int? local_disk_gb
  Int? num_preemptions

	command {
		set -euo pipefail

		if [ ${run_cosmo} == true ]; then
			echo "Running COSMO"
			
			if [ ${STANDALONE} == "false" ]; then
  			mkdir panoply_harmonize_output
  			tar -xf ${panoply_harmonize_tar} -C panoply_harmonize_output
  			tar_dir="$(pwd)/panoply_harmonize_output"
  			tar_name="$(ls panoply_harmonize_output | head -1)"
  
  			Rscript /prot/proteomics/Projects/PGDAC/src/cosmo/validate_harmonize_tar.R "$(pwd)/panoply_harmonize_output"
  
  			d1_file="$(cat data_file.txt)"
  			d2_file="$tar_dir/$tar_name/harmonized-data/rna-matrix.csv"
  			sample_file="$tar_dir/$tar_name/harmonized-data/sample-info.csv"
  	
    	else
        d1_file="${d1_file}"
  		  d2_file="${d2_file}"
  		  sample_file="${sample_file}"
  		fi

      yaml_file=${yaml_file}
			source /prot/proteomics/Projects/PGDAC/src/cosmo/panoply_run_cosmo.sh
			
			if [ ${STANDALONE} == false ]; then
		    cp $yaml_file $tar_dir/$tar_name/updated-master-parameter.yaml
			
			  cosmo_tar_dir="$tar_dir/$tar_name/cosmo-data"
			  mkdir $cosmo_tar_dir
			  cp -R $method1_out_folder $cosmo_tar_dir
        cp -R $method2_out_folder $cosmo_tar_dir
        cp -R $final_res_out_folder $cosmo_tar_dir
			  
			  home_dir=$(pwd)
			  cd $tar_dir
        tar -czvf "$home_dir/panoply_cosmo_output.tar" $tar_name
        cd $home_dir
        
      else
        tar -czvf panoply_cosmo_output.tar $out_dir
			
			fi

		else
			echo "COSMO not run"
			
			if [ ${STANDALONE} == false ]; then
			  mv ${panoply_harmonize_tar} panoply_cosmo_output.tar
			else
			  touch panoply_cosmo_output.tar
		  fi
			
		fi

		echo "Done"

	}

	runtime {
		docker: "broadcptacdev/panoply_cosmo:latest"
		memory: "${if defined(memory) then memory else '16'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
        cpu: "${if defined(cpu) then cpu else '6'}"
	}

	output {
		File cosmo_tar = "panoply_cosmo_output.tar"
		String d1_file_name = read_string("d1_file_name.txt")
		String d2_file_name = read_string("d2_file_name.txt")
	}
}
