#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#


task panoply_clumps_ptm_diffexp {
	# PTM GCT files; must include at least one
	File? pSTY_gct
	File? acK_gct
	File? ubK_gct

	String? gene_column
  	Float? fdr_assoc

	File groupsFile

	String output_prefix="results"
	File yaml_file

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail

		Rscript /prot/proteomics/Projects/PGDAC/src/clumps_diffexp.r ${"-p " + pSTY_gct} ${"-a " + acK_gct} ${"-u " + ubK_gct} \
		-g ${groupsFile} ${"-c " + gene_column} ${"-f " + fdr_assoc} -x ${output_prefix} ${"-y " + yaml_file} --libdir /prot/proteomics/Projects/PGDAC/src/ 
	}

	output {
		File results="${output_prefix}_full_mapped_sites_to_pdbs.tsv", # all var-sites with PDB mapping results
	}

	runtime {
		docker : "broadcptacdev/panoply_association:latest"
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

task panoply_clumps_ptm {
	# Mapping & Differential Expression Results
	File? de_file
	File? var_sites_file

	# ID Mapping
	File FASTA_ref_file			# file with FASTA reference sequences, to be blasted against UNIPROT
	String? accession_col		# rdesc column with protein accession id; must match IDs in FASTA_ref_file
	String? variable_sites_col	# rdesc column with variable sites (e.g. 'T527t')
	String? variable_sites_sep	# seperator for variable sites
	String? gene_column

	String annot
	String output_prefix="results"
	File yaml_file

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemtions
	
	command {
		set -euo pipefail


		mkdir clumpsptm_runs

		## for each subgroup in annotation

		for GROUP in (groups); do
			clumpsptm -i ${de_file} #--f phosphoproteome \
				-w Fold.Change --maps {var_sites_file} --pdbstore /pdbs \
				--grouping ${GROUP} --protein_id ${accession_col} \
				--threads 12 -v --subset positive --output_dir clumpsptm_runs/${GROUP}_pos_results
		done


		for GROUP in (groups); do
			clumpsptm -i /opt/input/pSTY_DE_analysis_K3_long_filt.tsv --f phosphoproteome \
				-w Fold.Change --maps /opt/input/mapped_sites_to_pdbs_noID.tsv --pdbstore /pdbs \
				--grouping ${GROUP} --protein_id id.description \
				--threads 12 -v --subset negative --output_dir clumpsptm_runs/${GROUP}_neg_results
		done


		## tar results 
	}

	output {
		File results="${output_prefix}_${annot}_clumpsResults.tar.gz"
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
	# PTM GCT files; must include at least one
	File? pSTY_gct
	File? acK_gct
	File? ubK_gct

	# ID Mapping
	File FASTA_ref_file			# file with FASTA reference sequences, to be blasted against UNIPROT
	String? accession_col		# rdesc column with protein accession id; must match IDs in FASTA_ref_file
	String? variable_sites_col	# rdesc column with variable sites (e.g. 'T527t')
	String? variable_sites_sep	# seperator for variable sites
	String? gene_column

	String output_prefix="results"
	File yaml_file


	# run ssGSEA on the geneset
	scatter (f in subset_files) {
		call panoply_clumps_ptm_diffexp {
			input:
				input_ds="${subset_bucket}/${f}",
				gene_set_database=panoply_cmap_input.genesets,
				yaml=yaml
		}
		call panoply_clumps_ptm {
			input:
				input_ds="${subset_bucket}/${f}",
				gene_set_database=panoply_cmap_input.genesets,
				yaml=yaml
		}
	}

    call panoply_clumps_ptm_mapping
}
