#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#

task panoply_clumps_ptm {
		File diff_exp_file
		File var_sites_file

		# List of Tarfiles in PDB directory, resolved by panoply_clumps_ptm_workflow.wdl
		Array[File]+ PDB_DIR

		String? accession_col						# id column (in var_sites_file / diff_exp_file) with protein accession id
		String? variable_sites_col					# column (in var_sites_file) with variable sites (e.g. 'T527t')

		Boolean? run_combined						# toggle for running all PTM sites combined

		String? weight_col							# column (in diff_exp_file) with weights to use for clumpsptm

		String output_prefix="results"
		File yaml_file

		Boolean DEBUG_MODE = false			# turn on debug mode, which limits the number of proteins mapped to 50 (randomly chosen)

		Int? memory
		Int? disk_space
		Int num_threads = 32 		# set default in inputs, rather than in runtime, so the argument can be used by clumps
		Int? num_preemptions

	command {
		set -euo pipefail

		# Unpack the PDB Archive
		echo "[`date +'%Y-%m-%d %T'`] INFO: Untarring PDB Archive"
		pdb_dir=pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/
		mkdir -p $pdb_dir # make PDB directory
		parallel -j ${num_threads} "tar -C $pdb_dir -xf" ::: ${sep=" " PDB_DIR} # untar each tar file
		echo "[`date +'%Y-%m-%d %T'`] INFO: Finished untarring PDB Archive"
		parallel -j ${num_threads} 'rm' ::: ${sep=' ' PDB_DIR} # remove tar-files to save space
		echo "[`date +'%Y-%m-%d %T'`] INFO: Finished removing PDB Tars"

		mkdir clumpsptm_runs

		python -u /prot/proteomics/Projects/PGDAC/src/clumps_ptm_wrapper.py -y ${yaml_file} \
			--input ${diff_exp_file} --maps ${var_sites_file} --pdbstore pdbs/ \
			${'--protein_id ' + accession_col} ${'--site_id ' + variable_sites_col} \
			${'--weight ' + weight_col} \
			${if defined(run_combined) then "--run_combined ${if select_first([run_combined]) then 'true' else 'false'}" else ""} \
			--threads ${num_threads} \
			$( [ ${DEBUG_MODE} = true ] && echo "-t" )

			# testing arguments
			# --protein_id ${accession_col} --weight ${weight_col} --threads ${num_threads} --run_combined ${run_combined} # for testing locally
			# group=1
			# --features phosphoproteome ubiquitylome acetylome --grouping $group --subset positive --output_dir "clumpsptm_runs/"$group"_pos_results"  # for testing clumpptm


		tar -czf ${output_prefix}_clumps_runs.tar -C clumpsptm_runs/ . # tar results
	}

	output {
		File results="${output_prefix}_clumps_runs.tar"
	}

	runtime {
		docker : "broadcptacdev/panoply_clumps_ptm:latest"
		memory : select_first ([memory, 32]) + "GB"
		disks : "local-disk  " + select_first ([disk_space, 100]) + " HDD"
		cpu : num_threads				# default set in inputs
		preemptible : select_first ([num_preemptions, 0])
	}

	meta {
		author : "C.M. Williams"
		email : "proteogenomics@broadinstitute.org"
	}

}

################################################
## workflow
workflow panoply_clumps_ptm_workflow {
	call panoply_clumps_ptm

}
