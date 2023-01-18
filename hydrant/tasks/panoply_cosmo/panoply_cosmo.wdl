workflow panoply_cosmo_workflow {
	String STANDALONE

	call panoply_cosmo {
		input:
			STANDALONE = STANDALONE
	}
	
	if(panoply_cosmo.run_cosmo_final) {
	  call panoply_cosmo_report {
	    input:
	      STANDALONE = STANDALONE,
	      cosmo_output_tar = panoply_cosmo.cosmo_tar,
	      d1_file_name = panoply_cosmo.d1_file_name,
        d2_file_name = panoply_cosmo.d2_file_name
	  }
	}
}

task panoply_cosmo {
	String STANDALONE
	File yaml_file
	Boolean run_cosmo = true
	String? sample_label

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
		
		echo "Managing parameters"

		Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
			--module cosmo \
			--master_yaml ${yaml_file} \
			--cosmo_run_cosmo ${run_cosmo} \
			${if defined(sample_label) then "--cosmo_sample_label " + "'${sample_label}'" else ""}

    yaml_file="final_output_params.yaml"
		R -e "cat(yaml::read_yaml('$yaml_file')[['cosmo.params']][['run_cosmo']], file = 'run_cosmo.txt', sep = '\n')"
		run_cosmo=$(cat 'run_cosmo.txt')

		if [ $run_cosmo == true ]; then
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
		  
		  touch d1_file_name.txt
			touch d2_file_name.txt
			
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
		Boolean run_cosmo_final = read_boolean("run_cosmo.txt")
		File yaml_file_final = "final_output_params.yaml"
		String d1_file_name = read_string("d1_file_name.txt")
		String d2_file_name = read_string("d2_file_name.txt")
	}
}

task panoply_cosmo_report {
  String STANDALONE
	File cosmo_output_tar
	String d1_file_name
	String d2_file_name

  Int? memory
  Int? local_disk_gb
  Int? num_preemptions

  command {
    set -euo pipefail
  
    tar_dir="$(pwd)/tar_output"
    mkdir $tar_dir
  	tar -xf ${cosmo_output_tar} -C $tar_dir
  	
  	if [ ${STANDALONE} == false ]; then
  	  cosmo_res_path="$tar_dir/$(ls $tar_dir | head -1)/cosmo-data/final_res_folder/cosmo_final_result.tsv"
  	  sample_corr_path="$tar_dir/$(ls $tar_dir | head -1)/cosmo-data/method1_folder/sample_correlation.csv"
  	else
  	  cosmo_res_path="$tar_dir/$(ls $tar_dir | head -1)/final_res_folder/cosmo_final_result.tsv"
  	  sample_corr_path="$tar_dir/$(ls $tar_dir | head -1)/method1_folder/sample_correlation.csv"
    fi

  	R -e \
  	  "rmarkdown::render('/prot/proteomics/Projects/PGDAC/src/cosmo/panoply_cosmo_report.Rmd', 
  	  params = list(final_result_path = '$cosmo_res_path', d1_file_name = '${d1_file_name}', d2_file_name = '${d2_file_name}', sample_corr_path = '$sample_corr_path'),
      output_dir = getwd())"
  }

  output {
      File cosmo_report_html = "panoply_cosmo_report.html"
  }

  runtime {
      docker : "broadcptacdev/panoply_cosmo:latest"
      memory: "${if defined(memory) then memory else '2'}GB"
      disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
      preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
  }

  meta {
      author : "Stephanie Vartany"
      email : "proteogenomics@broadinstitute.org"
  }
}