workflow panoply_omicsev_workflow {
    call OmicsEV_task
}

task OmicsEV_task {
    File panoply_harmonize_tar_file
    File yaml_file
	String? class_column_name

    Int? cpu
    Int? memory
    String? data_type
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
    ${default="protein" data_type} \
	x2.tsv
    }

    runtime {
    	docker: "broadcptacdev/panoply_omicsev:latest"
        memory: "${if defined(memory) then memory else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
        cpu: "${if defined(cpu) then cpu else '6'}"
    }
    output {
        File html_report = "final_evaluation_report.html"       
    }
}