#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#


task panoply_clumps_ptm {
	File diff_exp_file
	File var_sites_file

	File PDB_DIR				# PDB Directory tarfile

	String output_prefix="results"
	File yaml_file

	Int? memory
	Int? disk_space
	Int? num_threads=16 		# set default in inputs, rather than in runtime, so the argument can be used by clumps
	Int? num_preemtions
	
	command <<<
		set -euo pipefail
		
		# Unpack the PDB Archive
		mkdir -p /pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/
		pv ${PDB_DIR} | tar --use-compress-program=pigz -xf - -C /pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb

		mkdir clumpsptm_runs

		# drop first row (header) and get all unique values from ID column
		groups=`tail -n +2 ${diff_exp_file} | awk '{print $3}' | sort -u | sed 's/^"\(.*\)"$/\1/'` # NOTE: fails if we use command {}
		echo $groups # echo to ensure the file was read properly

		# run separately on positive and negative 
		for group in $groups; do
			set +e # allow errors, to prevent clumpsptm fails from ending script

			clumpsptm -i ${diff_exp_file} --features phosphoproteome ubiquitylome acetylome \
				-w logFC --maps ${var_sites_file} --pdbstore /pdbs \
				--grouping $group --protein_id id.description \
				--threads ${num_threads} -v --subset positive --output_dir "clumpsptm_runs/"$group"_pos_results" 2> err.txt
			
			# check if previous command had an error
			if [ $? -ne 0 ]; then
				# if the error message is NOT a "NO RESULTS FOUND" error, stop the script
				if [ -z `grep -l err.txt -e "ValueError: NO RESULTS FILES FOUND."` ]; then
					cat err.txt >> /dev/stderr
					exit 1
				else # otherwise just print a warning
					echo "WARNING: No results found for "$group"."
				fi
			fi
		done

		for group in $groups; do
			set +e # allow errors, to prevent clumpsptm fails from ending script

			clumpsptm -i ${diff_exp_file} --features phosphoproteome ubiquitylome acetylome \
				-w logFC --maps ${var_sites_file} --pdbstore /pdbs \
				--grouping $group --protein_id id.description \
				--threads ${num_threads} -v --subset negative --output_dir "clumpsptm_runs/"$group"_neg_results" 2> err.txt

			# check if previous command had an error
			if [ $? -ne 0 ]; then
				# if the error message is NOT a "NO RESULTS FOUND" error, stop the script
				if [ -z `grep -l err.txt -e "ValueError: NO RESULTS FILES FOUND."` ]; then
					cat err.txt >> /dev/stderr
					exit 1
				else # otherwise just print a warning
					echo "WARNING: No results found for "$group"."
				fi
			fi
		done

		set -euo pipefail # turn error catching back on

		tar -czf ${output_prefix}_clumps_runs.tar -C clumpsptm_runs/ . # tar results
	>>>

	output {
		File results="${output_prefix}_clumps_runs.tar" # all differential-expression files
	}

	runtime {
		docker : "broadcptacdev/panoply_clumps_ptm:latest"
		memory : select_first ([memory, 32]) + "GB"
		disks : "local-disk  " + select_first ([disk_space, 20]) + " HDD"
		cpu : num_threads				# default set in inputs
		preemptible : select_first ([num_preemtions, 0])
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
