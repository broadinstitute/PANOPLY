task panoply_cosmo_report {
	File cosmo_output

    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    command {
        set -euo pipefail

        unzip ${cosmo_output}

		R -e "rmarkdown::render('/prot/proteomics/Projects/PGDAC/src/panoply_cosmo_report.Rmd', 
		params = list(final_table_path = paste(getwd(), '/output/final_res_folder/cosmo_final_result.tsv', sep='')),
		output_dir = getwd())"
    }

    output {
        File cosmo_report_html = "panoply_cosmo_report.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_cosmo_report:latest"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Stephanie Vartany"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_cosmo_report_workflow {
    call panoply_cosmo_report
}