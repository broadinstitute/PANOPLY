#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#

task panoply_mice {
	File file

	Float? na_max
	Int? num_imputations
	Int? num_cores
	Int? random_seed

	String? output_prefix

	Int? memory
	Int? disk_space
	Int? num_preemptions
	
	command {
		set -euo pipefail
		
		Rscript /prot/proteomics/Projects/PGDAC/src/mice_imputation.R \
			${"--file_path " + file} \
			${"--na_max " + na_max} \
			${"--num_imputations " + num_imputations} \
			${"--num_cores " + num_cores} \
			${"--seed " + random_seed} \
			${"--output_prefix " + output_prefix} 
	}

	output {
	    File aggregate_imputation_gct = "${output_prefix}_mice_imputed.gct"
	    File aggregate_imputation_csv = "${output_prefix}_mice_imputed.csv"
	    File RDate_imputation_object = "${output_prefix}_mice_imputed.RData"
	}

    runtime {
	    docker: "broadcptacdev/panoply_mice:latest"
	    memory : select_first ([memory, 12]) + "GB"
	    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
	    cpu : select_first ([num_cores, 1]) + ""
	    preemptible : select_first ([num_preemptions, 0])
    }

	meta {
		author : "Stephanie Vartany"
		email : "proteogenomics@broadinstitute.org"
	}

}

################################################
## workflow
workflow panoply_mice_workflow {
    call panoply_mice
}
