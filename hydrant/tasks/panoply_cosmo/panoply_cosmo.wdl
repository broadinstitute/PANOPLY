#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

workflow panoply_cosmo_workflow {
	String STANDALONE
	File yaml_file
	File? panoply_harmonize_tar
	String label
	String? ome_type

	call panoply_cosmo {
		input:
			STANDALONE = STANDALONE,
			yaml_file = yaml_file,
			panoply_harmonize_tar = panoply_harmonize_tar,
			ome_type = ome_type,
	}
	
	if(panoply_cosmo.run_cosmo_final) {
	  call panoply_cosmo_report {
	    input:
	      cosmo_output_tar = panoply_cosmo.cosmo_tar,
	      d1_file_name = panoply_cosmo.d1_file_name,
        d2_file_name = panoply_cosmo.d2_file_name,
        label = label
	  }
	}
	
	output {
	  File cosmo_tar = panoply_cosmo.cosmo_tar
	  File? cosmo_report = panoply_cosmo_report.cosmo_report_html
	}
}

task panoply_cosmo {
	String STANDALONE
	File yaml_file
	Boolean? run_cosmo
	String? sample_label
	String? ome_type

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
			${if defined(run_cosmo) then "--cosmo_run_cosmo " + "'${run_cosmo}'" else ""} \
			${if defined(sample_label) then "--cosmo_sample_label " + "'${sample_label}'" else ""}

    yaml_file="final_output_params.yaml"
		R -e "cat(tolower(yaml::read_yaml('$yaml_file')[['cosmo.params']][['run_cosmo']]), file = 'run_cosmo.txt', sep = '\n')"
		run_cosmo=$(cat 'run_cosmo.txt')

		if [ $run_cosmo == true ]; then
			echo "Running COSMO"
			
			if [ ${STANDALONE} == "false" ]; then
  			tar -xf ${panoply_harmonize_tar}
  	  
    	  # get the root directory of the tar file
    	  tar -tf ${panoply_harmonize_tar} > all_files_in_tar.txt
    	  tar_dir=$(pwd)/$(basename $(head -n 1 all_files_in_tar.txt))
    	  
    	  Rscript /prot/proteomics/Projects/PGDAC/src/cosmo/validate_harmonize_tar.R $tar_dir ${ome_type}
  
  			d1_file="$tar_dir/harmonized-data/${ome_type}-matrix.csv"
  			d2_file="$tar_dir/harmonized-data/rna-matrix.csv"
  			sample_file="$tar_dir/harmonized-data/sample-info.csv"
  	
    	else
        d1_file="${d1_file}"
  		  d2_file="${d2_file}"
  		  sample_file="${sample_file}"
  		fi
      
      # run code for COSMO
      # the following variables are set in panoply_run_cosmo.sh:
        # $method1_out_folder
        # $method2_out_folder
        # $final_res_out_folder
        # $out_dir
        # $d1_file_name
        # $d2_file_name
			source /prot/proteomics/Projects/PGDAC/src/cosmo/panoply_run_cosmo.sh
			
			cp $yaml_file $out_dir/updated-master-parameter.yaml
			tar -czvf "panoply_cosmo_output.tar" $(basename $out_dir)

		else
			echo "COSMO not run"
			
			mkdir cosmo-data
			echo "COSMO not run" > cosmo-data/README.txt
			tar -czvf "panoply_cosmo_output.tar" cosmo-data
			
		  touch d1_file_name.txt
			touch d2_file_name.txt
			
		fi

		echo "Done"

	}

	runtime {
		docker: "broadcptacdev/panoply_cosmo:latest"
		memory: "${if defined(memory) then memory else '16'}GB"
    disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '32'} HDD"
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
	File cosmo_output_tar
	String d1_file_name
	String d2_file_name
	String label

  Int? memory
  Int? local_disk_gb
  Int? num_preemptions

  command {
    set -euo pipefail
  
  	tar -xf ${cosmo_output_tar}
  	
	  cosmo_res_path="$(pwd)/cosmo-data/final_res_folder/cosmo_final_result.tsv"
	  sample_corr_path="$(pwd)/cosmo-data/method1_folder/sample_correlation.csv"

  	R -e \
  	  "rmarkdown::render('/prot/proteomics/Projects/PGDAC/src/cosmo/panoply_cosmo_report.Rmd', 
  	  params = list(final_result_path = '$cosmo_res_path', d1_file_name = '${d1_file_name}', d2_file_name = '${d2_file_name}', sample_corr_path = '$sample_corr_path'),
      output_dir = getwd(),
      output_file = 'cosmo_${label}.html')"
  }

  output {
      File cosmo_report_html = "cosmo_" + label + ".html"
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