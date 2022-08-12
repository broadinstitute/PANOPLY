workflow panoply_omicsev_workflow {
    File data_file
    File sample_anno_file
    File? rna_file
    Int? cpu
    Int? memory
    String? data_type
    Int? local_disk_gb
    Int? num_preemptions

    call OmicsEV_task {
	input:
	    data_file=data_file, 
	    sample_anno_file=sample_anno_file, 
	    rna_file=rna_file,
	    data_type=data_type,
	    local_disk_gb=local_disk_gb,
	    num_preemptions= num_preemptions,
	    memory=memory,
	    cpu=cpu
    }

}

task OmicsEV_task {
    data_file=data_file, 
    sample_anno_file=sample_anno_file, 
    rna_file=rna_file,
    data_type=data_type,
    local_disk_gb=local_disk_gb,
    num_preemptions= num_preemptions,
    memory=memory,
    cpu=cpu

    command {
	set -euo pipefail
    
    	Rscript \
        /src/panoply_omicsev_preprocessing.R \
	${data_file} \
	${sample_anno_file} \
	${rna_file}
    
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