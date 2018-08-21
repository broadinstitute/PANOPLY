
task cmap_input {
  File tarball   # output from pgdac_cna_correlation
  String? cmap_group
  String? cmap_type
  String? cmap_log
  Int? cmap_permutations
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_cmapsetup-output.tar"
  String outGmtFile = "cmap-trans-genesets.gmt"
  
  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions
  
  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CMAPsetup -i ${tarball} -c ${codeDir} -o ${outFile} ${"-CMAPgroup " + cmap_group} ${"-CMAPtype " + cmap_type} ${"-CMAPnperm " + cmap_permutations} ${"-CMAPlog " + cmap_log}
  }
  
  output {
    File outputs = "${outFile}"
    File genesets = "${outGmtFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1.1"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}  



task cmap_ssgsea {
  # task adapted from pgdac_ssgsea; many inputs are set to specfic values for CMAP analysis
	File input_ds
	File gene_set_database
	String output_prefix = basename (input_ds, ".gctx")
	
	# other ssgsea options (below) are fixed for CMAP analysis
	String sample_norm_type = "rank"
	String correl_type = "rank"
	String statistic = "Kolmogorov-Smirnov"
	String output_score_type = "NES"

	Float weight = 0
	Int min_overlap = 5
	Int nperm = 0
	String global_fdr = "TRUE"
	String export_sigs = "FALSE"
	String ext_output = "FALSE"


  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions
	
	command {
		set -euo pipefail
		/home/pgdac/ssgsea-cli.R -i ${input_ds} -d ${gene_set_database} -o ${output_prefix} -n ${sample_norm_type} -w ${weight} -c ${correl_type} -t ${statistic} -s ${output_score_type} -p ${nperm} -m ${min_overlap} -g ${global_fdr} -e ${export_sigs} -x ${ext_output}
	}

	output {
		File scores="${output_prefix}-scores.gct"
	}

	runtime {
		docker : "broadcptac/pgdac_ssgsea:3"
    memory : select_first ([memory, 60]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, 16]) + ""
    preemptible : select_first ([num_preemptions, 2])
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}
}



task cmap_connectivity {
  File tarball
  Array[File] subset_scores
  String scores_dir = "cmap-subset-scores"
  String outFile = "pgdac_cmap-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    # setup output scores directory ...
    if [ ! -d ${scores_dir} ]; then 
      mkdir ${scores_dir} 
    fi
    # ... and copy subset scores
    mv ${sep=" " subset_scores} ${scores_dir}
    
    # combine shards/gather and run conectivity score calculations
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CMAPconn -i ${tarball} -o ${outFile} -CMAPscr ${scores_dir}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1.1"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}




workflow run_cmap_analysis {
  File CNAcorr_tarball
  File subsetListFile
  String subset_bucket
  Array[String] subset_files = read_lines ("${subsetListFile}")

  
  call cmap_input {
    input:
      tarball=CNAcorr_tarball
  }

  scatter (f in subset_files) {
    call cmap_ssgsea {
      input:
	      input_ds="${subset_bucket}/${f}",
	      gene_set_database=cmap_input.genesets
    }
  }
  
  call cmap_connectivity {
    input:
      tarball=cmap_input.outputs,
      subset_scores=cmap_ssgsea.scores
  }

  output {
    File output = cmap_connectivity.outputs
  }
}


