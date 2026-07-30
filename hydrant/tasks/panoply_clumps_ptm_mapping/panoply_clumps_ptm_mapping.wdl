#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
version 1.1

task panoply_clumps_ptm_mapping {
	input {
		# PTM GCT files; must include at least one
		File? pSTY_gct
		File? acK_gct
		File? ubK_gct

		# List of Tarfiles in PDB directory, resolved by panoply_clumps_ptm_workflow.wdl
		Array[File]+ PDB_DIR

		# ID Mapping
		File FASTA_ref_file				# file with FASTA reference sequences, to be blasted against UNIPROT
		String? FASTA_sep_type			# file with FASTA reference sequences, to be blasted against UNIPROT

		String? accession_col			# rdesc column with protein accession id; must match IDs in FASTA_ref_file
		String? gene_column

		String? variable_sites_col		# rdesc column with variable sites (e.g. 'T527t')
		String? variable_sites_sep		# seperator for variable sites
		Boolean? keep_multi_sites		# should multi-site PTMs be mapped to PDBs
		Boolean? filter_duplicate_sites	# should multi-site PTMs that were also observed as single-sites be filtered out

		File UNIPROT_SWISSPROT			# file with UNIPROT sequences, to be blasted to
		File SIFTS_DB					# SIFTS database with mapping between UNIPROT and PDB IDs



		String output_prefix="results"
		File yaml_file

		Boolean DEBUG_MODE = false			# turn on debug mode, which limits the number of proteins mapped to 50 (randomly chosen)

		Int? memory
		Int? disk_space
		Int num_threads = 32 		# set default in inputs, rather than in runtime, so the argument can be used by clumps
		Int? num_preemptions
	}

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

		# Run Mapping Scripts 
		python /prot/proteomics/Projects/PGDAC/src/clumps_ptm_mapping.py --PDB_DIR pdbs \
		${'--phosphoproteome_gct ' + pSTY_gct} ${'--acetylome_gct ' + acK_gct} ${'--ubiquitylome_gct ' + ubK_gct} \
		--FASTA_ref_file ${FASTA_ref_file} ${'--FASTA_sep_type ' + FASTA_sep_type}  \
		${'--accession_col ' + accession_col} ${'--gene_column ' + gene_column} \
		${'--variable_sites_col ' + variable_sites_col} ${'--variable_sites_sep ' + '"' + variable_sites_sep + '"'} \
		${if defined(keep_multi_sites) then "--keep_multi_sites " + (if select_first([keep_multi_sites]) then "true" else "false") else ""} \
		${if defined(filter_duplicate_sites) then "--filter_duplicate_sites " + (if select_first([filter_duplicate_sites]) then "true" else "false") else ""} \
		--UNIPROT_SWISSPROT ${UNIPROT_SWISSPROT} --SIFTS_DB ${SIFTS_DB} \
		--output_prefix ${output_prefix} --yaml ${yaml_file} --num_threads ${num_threads} \
		$( [ ${DEBUG_MODE} = true ] && echo "--DEBUG_MODE" )
	}

	output {
		File full_results = "output_files/${output_prefix}_full_mapped_sites_to_pdbs.tsv" # all var-sites with PDB mapping results
		File filt_results = "output_files/${output_prefix}_mapped_sites_to_pdbs.tsv" # var-sites, filtered to valid PDB mappings
		File var_sites_file = "output_files/${output_prefix}_var_sites_combined.tsv" # file with all variable sites
		File mapping_params = "output_files/params.yaml" # parameters file
	}

	runtime {
		docker : "broadcptacdev/panoply_clumps_ptm_mapping:latest"
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
workflow panoply_clumps_ptm_mapping_workflow {
    call panoply_clumps_ptm_mapping

}