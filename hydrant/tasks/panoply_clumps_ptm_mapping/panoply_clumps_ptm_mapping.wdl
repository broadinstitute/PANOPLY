#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#

task panoply_clumps_ptm {
	# PTM GCT files; must include at least one
	File? pSTY_gct
	File? acK_gct
	File? ubK_gct

	File PDB_DIR				# PDB Directory tarfile

	# ID Mapping
	File FASTA_ref_file			# file with FASTA reference sequences, to be blasted against UNIPROT
	String? accession_col		# rdesc column with protein accession id; must match IDs in FASTA_ref_file
	String? variable_sites_col	# rdesc column with variable sites (e.g. 'T527t')
	String? variable_sites_sep	# seperator for variable sites
	String? gene_column


	String output_prefix="results"
	File yaml_file

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail

		# Unpack the PDB Archive
		mkdir -p /pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/
		pv ${PDB_DIR} | tar --use-compress-program=pigz -xf - -C /pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb
		#pv /pdb_folder/pdbs.tar.gz | tar -xzf - -C /pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb
		#pv /pdb_folder/pdbs.tar.zst | tar --zstd -xf - -C /pdbs/ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb

		## Run Mapping Scriptls 
		#python /prot/proteomics/Projects/PGDAC/src/clumps_ptm_mapping.py -p /opt/input/phosphoproteome-subset.gct -f /opt/input/Ensembl.human.hg19.clean3nr.602contams_20230913.fasta -i id.description -o ODG_v3 -y /opt/input/master-parameters.yaml -n 12 -d

		# Run Mapping Scriptls 
		python /prot/proteomics/Projects/PGDAC/src/clumps_ptm_mapping.py ${'--phosphoproteome_gct ' + pSTY_gct} ${'--acetylome_gct ' + acK_gct} ${'--ubiquitylome_gct ' + ubK_gct} '--PDB_DIR ' + $PDB_DIR '--FASTA_ref_file ' + FASTA_ref_file ${'--accession_col ' + accession_col} ${'--variable_sites_col ' + variable_sites_col} ${'--variable_sites_sep ' + '"' + variable_sites_sep + '"'} ${'--gene_column ' + gene_column} ${'--output_prefix ' + output_prefix} "--yaml " + yaml_file "--num_threads " + num_threads --libdir /prot/proteomics/Projects/PGDAC/src/
	}

	output {
		File full_results="${output_prefix}_full_mapped_sites_to_pdbs.tsv", # all var-sites with PDB mapping results
		File filt_results="output_files/${output_prefix}_mapped_sites_to_pdbs.tsv" # var-sites, filtered to valid PDB mappings
		File var_sites_file="output_files/${output_prefix}_var_sites_combined" # file with all variable sites
	}

	runtime {
		docker : "broadcptacdev/panoply_clumps_ptm_mapping:latest"
		memory : select_first ([memory, 32]) + "GB"
		disks : "local-disk  " + select_first ([disk_space, 20]) + " HDD"
		cpu : select_first ([num_threads, 32]) + ""
		preemptible : select_first ([num_preemtions, 0])
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
