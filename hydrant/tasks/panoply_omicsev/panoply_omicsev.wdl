workflow myWorkflow {
    File data_file
    File sample_anno_file
    File? rna_file
    Int? cpu
    Int? memory
    String? data_type
    Int? local_disk_gb
    Int? num_preemptions
    

    call OmicsEV_preprocessing {
	input:
	    data_file=data_file, 
	    sample_anno_file=sample_anno_file, 
	    rna_file=rna_file,
	    local_disk_gb=local_disk_gb,
	    num_preemptions= num_preemptions
    }

    call OmicsEV_task {
	input: 
	    data_zip=OmicsEV_preprocessing.data_zip, 
	    sample_list=OmicsEV_preprocessing.sample_list, 
	    x2=OmicsEV_preprocessing.x2, 
	    cpu=cpu, 
	    data_type=data_type, 
            code_file=OmicsEV_preprocessing.code_file,
            memory=memory,
	    local_disk_gb=local_disk_gb,
	    num_preemptions= num_preemptions
    }
}

task OmicsEV_preprocessing {
    File data_file
    File sample_anno_file
    File? rna_file
    Int? local_disk_gb
    Int? num_preemptions

    command {
	set -euo pipefail

	Rscript \
        /prot/proteomics/Projects/PGDAC/src/panoply_omicsev_preprocessing.R \
	${data_file} \
	${sample_anno_file} \
	${rna_file}
    }

    runtime {
	docker: "broadcptacdev/panoply_omicsev:latest"
	memory: "32GB"
	disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    output {
	File data_zip = "dataset.zip"
	File sample_list = "sample_list.tsv"
	File x2 = "x2.tsv"
	File code_file = "/prot/proteomics/Projects/PGDAC/src/panoply_run_OmicsEV.R"
    }
}

task OmicsEV_task {
    File data_zip
    File sample_list
    File x2
    Int cpu
    Int memory
    String? data_type
    File code_file
    Int? local_disk_gb
    Int? num_preemptions

    command {
	set -euo pipefail
    
    	unzip ${data_zip}
    
        Rscript ${code_file} \
        ${data_zip} \
        ${sample_list} \
        ${x2} \
        ${default=6 cpu} \
        ${default="protein" data_type}
    }

    runtime {
    	docker: "proteomics/omicsev:latest"
        memory: "${if defined(memory) then memory else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
        cpu: cpu
    }
    output {
        File html_report = "final_evaluation_report.html"       
    }
}