#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cmap_connectivity {
  File tarball
  File yaml
  String? cmap_grp
  String? cmap_typ
  Int? permutations
  Array[File] subset_scores
  Array[File]? permutation_scores
  String scores_dir = "cmap-subset-scores"
  String permutation_dir = "cmap-permutation-scores"
  String outFile = "panoply_cmap-output.tar"

  Float? cna_threshold
  Int? cna_effects_threshold
  Int? min_sigevents
  Int? max_sigevents
  Int? top_N
  Float? fdr_pvalue
  String? log_transform
  String? must_include_genes
  String? cis_fdr
  String? legacy_score
  Int? rankpt_n
  Int? mean_rankpt_threshold
  Float? cmap_fdr
  String? alpha


  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  String cmap_group = "${if defined (cmap_grp) then cmap_grp else 'all'}"
  String cmap_type = "${if defined (cmap_typ) then cmap_typ else 'pome'}"

  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module cmap_analysis \
    --master_yaml ${yaml} \
    ${"--cna_threshold " + cna_threshold} \
    ${"--cna_effects_threshold " + cna_effects_threshold} \
    ${"--min_sigevents " + min_sigevents} \
    ${"--max_sigevents " + max_sigevents} \
    ${"--top_N " + top_N} \
    ${"--fdr_pvalue " + fdr_pvalue} \
    ${"--log_transform " + log_transform} \
    ${"--must_include_genes " + must_include_genes} \
    ${"--cis_fdr " + cis_fdr} \
    ${"--legacy_score " + legacy_score} \
    ${"--rankpt_n " + rankpt_n} \
    ${"--mean_rankpt_threshold " + mean_rankpt_threshold} \
    ${"--cmap_fdr " + cmap_fdr} \
    ${"--alpha " + alpha}

    # setup output scores directory ...
    if [ ! -d ${scores_dir} ]; then
      mkdir ${scores_dir}
    fi
    # ... and copy subset scores
    mv ${sep=" " subset_scores} ${scores_dir}

    # same for perumations scores ...
    if [ ${permutations} -gt 0 ]; then
      if [ ! -d ${permutation_dir} ]; then
        mkdir ${permutation_dir}
      fi
      mv ${sep=" " permutation_scores} ${permutation_dir}
    fi

    # combine shards/gather and run conectivity score calculations
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CMAPconn -i ${tarball} -o ${outFile} -CMAPscr ${scores_dir} -CMAPnperm ${default="0" permutations} -CMAPpmt ${permutation_dir} -CMAPcfg "cmap-config-custom.r" ${"-CMAPgroup " + cmap_group} ${"-CMAPtype " + cmap_type}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptacdev/panoply_cmap_analysis:latest"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, permutations+1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task panoply_cmap_input {
  File tarball   # output from panoply_cna_correlation
  File yaml
  String? cmap_grp
  String? cmap_typ
  Int? cmap_permutations
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_cmapsetup-output.tar"
  String outGmtFile = "cmap-trans-genesets.gmt"

  Float? cna_threshold
  Int? cna_effects_threshold
  Int? min_sigevents
  Int? max_sigevents
  Int? top_N
  Float? fdr_pvalue
  String? log_transform
  String? must_include_genes
  String? cis_fdr
  String? legacy_score
  Int? rankpt_n
  Int? mean_rankpt_threshold
  Float? cmap_fdr
  String? alpha

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  String cmap_group = "${if defined (cmap_grp) then cmap_grp else 'all'}"
  String cmap_type = "${if defined (cmap_typ) then cmap_typ else 'pome'}"

  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module cmap_analysis \
    --master_yaml ${yaml} \
    ${"--cna_threshold " + cna_threshold} \
    ${"--cna_effects_threshold " + cna_effects_threshold} \
    ${"--min_sigevents " + min_sigevents} \
    ${"--max_sigevents " + max_sigevents} \
    ${"--top_N " + top_N} \
    ${"--fdr_pvalue " + fdr_pvalue} \
    ${"--log_transform " + log_transform} \
    ${"--must_include_genes " + must_include_genes} \
    ${"--cis_fdr " + cis_fdr} \
    ${"--legacy_score " + legacy_score} \
    ${"--rankpt_n " + rankpt_n} \
    ${"--mean_rankpt_threshold " + mean_rankpt_threshold} \
    ${"--cmap_fdr " + cmap_fdr} \
    ${"--alpha " + alpha}
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CMAPsetup -i ${tarball} -c ${codeDir} -o ${outFile} ${"-CMAPgroup " + cmap_group} ${"-CMAPtype " + cmap_type} ${"-CMAPnperm " + cmap_permutations} -CMAPcfg "cmap-config-custom.r"
  }

  output {
    File outputs = "${outFile}"
    File genesets = "${outGmtFile}"
    Array[File] permuted_genesets = glob ("*-permuted-genes-*.gmt")
  }

  runtime {
    docker : "broadcptacdev/panoply_cmap_analysis:latest"
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


task panoply_cmap_annotate {
  File tarball                  # output from pgdac_cmap_connectivity
  File cmap_data_file           # CMAP level 5 geneKD data (gctx)
  File? cmap_enrichment_groups   # groups file (ala experiment design file)
  File yaml
  String? cmap_grp
  String? cmap_typ
  String outFile = "panoply_cmap-annotate-output.tar"

  Float? cna_threshold
  Int? cna_effects_threshold
  Int? min_sigevents
  Int? max_sigevents
  Int? top_N
  Float? fdr_pvalue
  String? log_transform
  String? must_include_genes
  String? cis_fdr
  String? legacy_score
  Int? rankpt_n
  Int? mean_rankpt_threshold
  Float? cmap_fdr
  String? alpha

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  String cmap_group = "${if defined (cmap_grp) then cmap_grp else 'all'}"
  String cmap_type = "${if defined (cmap_typ) then cmap_typ else 'pome'}"


  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module cmap_analysis \
    --master_yaml ${yaml} \
    ${"--cna_threshold " + cna_threshold} \
    ${"--cna_effects_threshold " + cna_effects_threshold} \
    ${"--min_sigevents " + min_sigevents} \
    ${"--max_sigevents " + max_sigevents} \
    ${"--top_N " + top_N} \
    ${"--fdr_pvalue " + fdr_pvalue} \
    ${"--log_transform " + log_transform} \
    ${"--must_include_genes " + must_include_genes} \
    ${"--cis_fdr " + cis_fdr} \
    ${"--legacy_score " + legacy_score} \
    ${"--rankpt_n " + rankpt_n} \
    ${"--mean_rankpt_threshold " + mean_rankpt_threshold} \
    ${"--cmap_fdr " + cmap_fdr} \
    ${"--alpha " + alpha}
    Rscript /prot/proteomics/Projects/PGDAC/src/cmap-annotate.R  ${tarball} ${cmap_data_file} ${cmap_group} ${cmap_type} ${cmap_enrichment_groups} ${outFile} "cmap-config-custom.r"
  }

  output {
    File outputs = "${outFile}"
    File gsea_input = "${cmap_group}-cmap-${cmap_type}-gsea-input.gct"
  }

  runtime {
    docker : "broadcptacdev/panoply_cmap_annotate:latest"
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



task panoply_cmap_ssgsea {
  # task adapted from panoply_ssgsea; many inputs are set to specfic values for CMAP analysis
	File input_ds
	File gene_set_database
	Int? permutation_num
	String output_prefix = "${basename (input_ds, '.gctx')}" + "${if defined (permutation_num) then '-'+permutation_num else ''}"

	# other ssgsea options (below) are fixed for CMAP analysis
	String sample_norm_type = "rank"
	String correl_type = "rank"
	String statistic = "Kolmogorov-Smirnov"
	String output_score_type = "NES"

	Float weight = 0.0
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
		docker : "broadcptacdev/panoply_ssgsea:latest"
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



task panoply_cmap_annotate_ssgsea {
  # task adapted from pgdac_ssgsea; many inputs are set to specific values for CMAP annotation
	File input_ds
	File gene_set_database
  	String outFile="panoply_cmap_annotate-ssgsea.tar"
  	String output_prefix = "${basename (input_ds, '.gct')}"

	# other ssgsea options (below) are fixed for CMAP analysis
	String sample_norm_type = "rank"
	String correl_type = "rank"
	String statistic = "Kolmogorov-Smirnov"
	String output_score_type = "NES"

	Float weight = 0.0
	Int min_overlap = 5
	Int nperm = 1000
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

    # assemble results into tar file
    find * -regextype posix-extended -regex "^signature_gct/.*.gct$|^${output_prefix}.*.gct$|^.*.log.txt$|^.*parameters.txt$" -print0 | tar -cvf ${outFile} --null -T -
	}

	output {
		File outputs="${outFile}"
	}

	runtime {
		docker : "broadcptacdev/panoply_ssgsea:latest"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, 16]) + ""
    preemptible : select_first ([num_preemptions, 2])
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}
}


workflow run_cmap_analysis {
  File CNAcorr_tarball
  File subsetListFile
  File cmap_level5_data
  File annotation_pathway_db
  String subset_bucket
  Int? n_permutations
  Array[String] subset_files = read_lines ("${subsetListFile}")
  String? group
  String? data_type
  File yaml
  File? cmap_enrichment_groups

  Float? cna_threshold
  Int? cna_effects_threshold
  Int? min_sigevents
  Int? max_sigevents
  Int? top_N
  Float? fdr_pvalue
  String? log_transform
  String? must_include_genes
  String? cis_fdr
  String? legacy_score
  Int? rankpt_n
  Int? mean_rankpt_threshold
  Float? cmap_fdr
  String? alpha


  call panoply_cmap_input {
    input:
      tarball=CNAcorr_tarball,
      cmap_permutations=n_permutations,
      cmap_grp=group,
      cmap_typ=data_type,
      yaml=yaml,
      cna_threshold=cna_threshold,
      cna_effects_threshold=cna_effects_threshold,
      min_sigevents=min_sigevents,
      max_sigevents=max_sigevents,
      top_N=top_N,
      fdr_pvalue=fdr_pvalue,
      log_transform=log_transform,
      must_include_genes=must_include_genes,
      cis_fdr=cis_fdr,
      legacy_score=legacy_score,
      rankpt_n=rankpt_n,
      mean_rankpt_threshold=mean_rankpt_threshold,
      cmap_fdr=cmap_fdr,
      alpha=alpha
  }

  # run ssGSEA on the geneset
  scatter (f in subset_files) {
    call panoply_cmap_ssgsea {
      input:
	      input_ds="${subset_bucket}/${f}",
	      gene_set_database=panoply_cmap_input.genesets
    }
  }

  # run ssGSEA on the permuted genesets (if any)
  if (n_permutations > 0) {
    Array[Pair[String,Int]] fxp = cross (subset_files, range (n_permutations))
    scatter (x in fxp) {
      call panoply_cmap_ssgsea as permutation {
        input:
          input_ds="${subset_bucket}/${x.left}",
	        gene_set_database=panoply_cmap_input.permuted_genesets[x.right],
	        permutation_num=x.right
      }
    }
  }

  call panoply_cmap_connectivity {
    input:
      tarball=panoply_cmap_input.outputs,
      subset_scores=panoply_cmap_ssgsea.scores,
      permutations=n_permutations,
      permutation_scores=permutation.scores,
      cmap_grp=group,
      cmap_typ=data_type,
      yaml=yaml,
      cna_threshold=cna_threshold,
      cna_effects_threshold=cna_effects_threshold,
      min_sigevents=min_sigevents,
      max_sigevents=max_sigevents,
      top_N=top_N,
      fdr_pvalue=fdr_pvalue,
      log_transform=log_transform,
      must_include_genes=must_include_genes,
      cis_fdr=cis_fdr,
      legacy_score=legacy_score,
      rankpt_n=rankpt_n,
      mean_rankpt_threshold=mean_rankpt_threshold,
      cmap_fdr=cmap_fdr,
      alpha=alpha
  }

  call panoply_cmap_annotate {
    input:
      tarball=panoply_cmap_connectivity.outputs,
      cmap_data_file=cmap_level5_data,
      cmap_grp=group,
      cmap_typ=data_type,
      yaml=yaml,
      cmap_enrichment_groups=cmap_enrichment_groups,
      cna_threshold=cna_threshold,
      cna_effects_threshold=cna_effects_threshold,
      min_sigevents=min_sigevents,
      max_sigevents=max_sigevents,
      top_N=top_N,
      fdr_pvalue=fdr_pvalue,
      log_transform=log_transform,
      must_include_genes=must_include_genes,
      cis_fdr=cis_fdr,
      legacy_score=legacy_score,
      rankpt_n=rankpt_n,
      mean_rankpt_threshold=mean_rankpt_threshold,
      cmap_fdr=cmap_fdr,
      alpha=alpha
  }

  call panoply_cmap_annotate_ssgsea {
    input:
      input_ds=panoply_cmap_annotate.gsea_input,
  	  gene_set_database=annotation_pathway_db
  }

  output {
    File outputs = panoply_cmap_annotate.outputs
    File ssgseaOutput = panoply_cmap_annotate_ssgsea.outputs
  }
}
