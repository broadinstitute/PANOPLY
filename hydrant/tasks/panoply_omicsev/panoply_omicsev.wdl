workflow panoply_omicsev_workflow {
	Boolean STANDALONE
    File yaml_file
    String? class_column_name
    
    Int? cpu
    Int? memory
    Int? local_disk_gb
    Int? num_preemptions

    if(STANDALONE) {
    	call OmicsEV_task_standalone {
        	input: 
            	yaml_file = yaml_file,
            	class_column_name = class_column_name,
                cpu = cpu,
                memory = memory,
                local_disk_gb = local_disk_gb,
                num_preemptions = num_preemptions
        }
    }
    
    if(!STANDALONE) {
    	call OmicsEV_task {
        	input: 
            	yaml_file = yaml_file,
            	class_column_name = class_column_name,
                cpu = cpu,
                memory = memory,
                local_disk_gb = local_disk_gb,
                num_preemptions = num_preemptions
        }
    }
}

task OmicsEV_task_standalone {
    Array[File]? data_files
    File? sample_anno_file
    File? rna_file
    File yaml_file
    String? class_column_name
    String? data_type
    
    Int? cpu
    Int? memory
    Int? local_disk_gb
    Int? num_preemptions

    command {
	set -euo pipefail
    
    Rscript \
    /src/panoply_omicsev_preprocessing_standalone.R \
	${sep=',' data_files} \
	${sample_anno_file} \
	${default = 'no_rna' rna_file} \
	${default = 'Type' class_column_name} \
	${yaml_file}
    
    Rscript \
	/src/panoply_run_OmicsEV.R \
    datasets \
    sample_list.tsv \
    ${default=6 cpu} \
    ${default="protein" data_type} \
	x2.tsv
    }

    runtime {
    	docker: "broadcptacdev/panoply_omicsev:latest"
        memory: "${if defined(memory) then memory else '96'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
        cpu: "${if defined(cpu) then cpu else '6'}"
    }
    output {
        File html_report = "final_evaluation_report.html"       
    }
    parameter_meta {
    	data_files: "Required for STANDALONE. Array of gct files to be used for OmicsEV. Contain same samples. Can be a single data file."
    	sample_anno_file: "Required for STANDALONE. Sample annotation csv file. Must have a column header 'Type' or the input for class_column_name."
        rna_file: "Optional. RNA gct file for correlation report."
        yaml_file: "Required. Yaml file in PANOPLY format."
        class_column_name: "Optional. A column header from sample_anno_file. Default is 'Type'."
        data_type: "Optional. The data type of the input file(s), either 'protein' or 'gene'. Default is 'protein'."
    }
}

task OmicsEV_task {
    File? panoply_harmonize_tar_file
    File yaml_file
    String? class_column_name
    
    Int? cpu
    Int? memory

    Int? local_disk_gb
    Int? num_preemptions

    command {
	set -euo pipefail
    
    Rscript \
    /src/panoply_omicsev_preprocessing.R \
	${panoply_harmonize_tar_file} \
	${yaml_file} \
	${default="Type" class_column_name}
    
    Rscript \
	/src/panoply_run_OmicsEV.R \
    dataset \
    sample_list.tsv \
    ${default=6 cpu} \
    protein \
	x2.tsv
    }

    runtime {
    	docker: "broadcptacdev/panoply_omicsev:latest"
        memory: "${if defined(memory) then memory else '96'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
        cpu: "${if defined(cpu) then cpu else '6'}"
    }
    output {
        File html_report = "final_evaluation_report.html"       
    }
    parameter_meta {
    	panoply_harmonize_tar_file: "Required if not STANDALONE. Tar file output from panoply_harmonize."
    	yaml_file: "Required. Yaml file in PANOPLY format."
        class_column_name: "Optional. A column header from sample_anno_file if the desired phemotype class is not 'Type'."
    }
}