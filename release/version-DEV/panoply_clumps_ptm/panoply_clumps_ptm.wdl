#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#

task panoply_clumps_ptm {
	File diff_exp_file
	File var_sites_file

	# Google-Cloud Bucket with PDB Directory split into tarfiles
	String PDB_ref_bucket

	# List of Tarfiles in PDB directory
	Array[File]+ PDB_DIR = [ PDB_ref_bucket + "/pdbs_0.tar", PDB_ref_bucket + "/pdbs_1.tar", PDB_ref_bucket + "/pdbs_2.tar", PDB_ref_bucket + "/pdbs_3.tar", PDB_ref_bucket + "/pdbs_4.tar", PDB_ref_bucket + "/pdbs_5.tar", PDB_ref_bucket + "/pdbs_6.tar", PDB_ref_bucket + "/pdbs_7.tar", PDB_ref_bucket + "/pdbs_8.tar", PDB_ref_bucket + "/pdbs_9.tar", PDB_ref_bucket + "/pdbs_a.tar", PDB_ref_bucket + "/pdbs_b.tar", PDB_ref_bucket + "/pdbs_c.tar", PDB_ref_bucket + "/pdbs_d.tar", PDB_ref_bucket + "/pdbs_e.tar", PDB_ref_bucket + "/pdbs_f.tar", PDB_ref_bucket + "/pdbs_g.tar", PDB_ref_bucket + "/pdbs_h.tar", PDB_ref_bucket + "/pdbs_i.tar", PDB_ref_bucket + "/pdbs_j.tar", PDB_ref_bucket + "/pdbs_k.tar", PDB_ref_bucket + "/pdbs_l.tar", PDB_ref_bucket + "/pdbs_m.tar", PDB_ref_bucket + "/pdbs_n.tar", PDB_ref_bucket + "/pdbs_o.tar", PDB_ref_bucket + "/pdbs_p.tar", PDB_ref_bucket + "/pdbs_q.tar", PDB_ref_bucket + "/pdbs_r.tar", PDB_ref_bucket + "/pdbs_s.tar", PDB_ref_bucket + "/pdbs_t.tar", PDB_ref_bucket + "/pdbs_u.tar", PDB_ref_bucket + "/pdbs_v.tar", PDB_ref_bucket + "/pdbs_w.tar", PDB_ref_bucket + "/pdbs_x.tar", PDB_ref_bucket + "/pdbs_y.tar", PDB_ref_bucket + "/pdbs_z.tar" ]

	String? accession_col						# id column (in var_sites_file / diff_exp_file) with protein accession id
	String? variable_sites_col					# column (in var_sites_file) with variable sites (e.g. 'T527t')

	Boolean? run_combined						# toggle for running all PTM sites combined

	String? weight_col							# column (in diff_exp_file) with weights to use for clumpsptm

	String output_prefix="results"
	File yaml_file

	Boolean? DEBUG_MODE=false			# turn on debug mode, which limits the number of proteins mapped to 50 (randomly chosen)

	Int? memory
	Int? disk_space
	Int? num_threads=32 		# set default in inputs, rather than in runtime, so the argument can be used by clumps
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
			${true="--run_combined true" false="--run_combined false" run_combined} \
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
		docker : "broadcptacdev/panoply_clumps_ptm:DEV"
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
