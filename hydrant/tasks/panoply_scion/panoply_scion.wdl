version 1.1

## Runs SCION network inference, sharding permutation-based FDR thresholding
## across one scatter shard per permutation instead of one long serial task.
## Backed by the `SCION` R package (github.com/nmclark2/SCION) rather than the
## loose scripts this task previously vendored -- see that repo's inst/wdl/
## for the same workflow packaged for standalone (non-PANOPLY) use, sharing
## this same Docker image.
##
## clustering_method defaults to "kmeans" with no clustering_data_file, matching
## the old task's always-cluster-via-kmeans default: cluster_genes() auto-derives
## the clustering matrix from the target+regulator data (regulator data takes
## precedence for genes that are both) when none is supplied.
##
## pome_gct_file's row names are rebuilt from its rdesc metadata (see
## scion_prepare_gct.R) rather than used as raw GCT row IDs: gene_id_col names
## the rdesc column holding the gene symbol (matching the gene_id_col
## convention used elsewhere in PANOPLY, e.g. panoply_harmonize), and for PTM
## omes (ome = phosphoproteome/acetylome/ubiquitylome) the PTM site is also
## extracted from the raw row ID and appended via ptm_sep, per the
## SpectrumMill (ptm_type = "SM") or FragPipe (ptm_type = "FP") ID convention.

workflow panoply_scion_workflow {
  input {
    String ome
    File pome_gct_file
    File mrna_gct_file
    File TF_file

    String gene_id_col = "geneSymbol"
    String ptm_type = "SM"
    String clustering_method = "kmeans"
    Float weight_threshold = 0
    Boolean normalize = false
    Int num_cores = 1
    String ptm_sep = "_"
    Int seed = 2023
    Int nb_trees = 10000

    Int n_permutations = 100
    String permute_dim = "col"
    Int base_seed = 0
    Float target_fdr = 0.05

    Int memory = 32
    Int disk_space = 50
    Int num_preemptions = 0
  }

  call panoply_scion_run_real {
    input:
      ome = ome,
      mrna_gct_file = mrna_gct_file,
      pome_gct_file = pome_gct_file,
      TF_file = TF_file,
      gene_id_col = gene_id_col,
      ptm_type = ptm_type,
      clustering_method = clustering_method,
      weight_threshold = weight_threshold,
      normalize = normalize,
      num_cores = num_cores,
      ptm_sep = ptm_sep,
      seed = seed,
      nb_trees = nb_trees,
      memory = memory,
      disk_space = disk_space,
      num_preemptions = num_preemptions
  }

  scatter (i in range(n_permutations)) {
    call panoply_scion_run_permutation {
      input:
        target_rds = panoply_scion_run_real.target_rds,
        reg_rds = panoply_scion_run_real.reg_rds,
        cluster_assignment_rds = panoply_scion_run_real.cluster_assignment_rds,
        params_rds = panoply_scion_run_real.params_rds,
        index = i + 1,
        base_seed = base_seed,
        permute_dim = permute_dim,
        num_preemptions = num_preemptions
    }
  }

  call panoply_scion_aggregate_fdr {
    input:
      network_rds = panoply_scion_run_real.network_rds,
      permutation_files = panoply_scion_run_permutation.permutation_rds,
      target_fdr = target_fdr,
      num_preemptions = num_preemptions
  }

  output {
    File real_network_tsv = panoply_scion_run_real.network_tsv
    File thresholded_network_tsv = panoply_scion_aggregate_fdr.thresholded_network_tsv
    File fdr_curve_png = panoply_scion_aggregate_fdr.fdr_curve_png
    File weight_comparison_png = panoply_scion_aggregate_fdr.weight_comparison_png
    File? network_plot_png = panoply_scion_aggregate_fdr.network_plot_png
    File fdr_result_rds = panoply_scion_aggregate_fdr.fdr_result_rds
  }

  meta {
    author : "Natalie Clark"
    email : "nclark@broadinstitute.org"
  }
}

task panoply_scion_run_real {
  input {
    String ome
    File mrna_gct_file
    File pome_gct_file
    File TF_file
    String gene_id_col
    String ptm_type
    String clustering_method
    Float weight_threshold
    Boolean normalize
    Int num_cores
    String ptm_sep
    Int seed
    Int nb_trees
    Int memory
    Int disk_space
    Int num_preemptions
  }

  command <<<
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/scion_run_real.R \
      --target_data_file ~{mrna_gct_file} \
      --reg_data_file ~{pome_gct_file} \
      --ome ~{ome} \
      --gene_id_col '~{gene_id_col}' \
      --ptm_type ~{ptm_type} \
      --reg_genes_file ~{TF_file} --gene_list_header FALSE \
      --clustering_method ~{clustering_method} \
      --weightthreshold ~{weight_threshold} \
      --normalize ~{normalize} \
      --num_cores ~{num_cores} \
      --ptm_sep '~{ptm_sep}' \
      --seed ~{seed} \
      --nb_trees ~{nb_trees} \
      --out_dir out

    cp out/network.tsv "out/~{ome}-SCION-network-full.tsv"
  >>>

  output {
    File network_tsv = "out/~{ome}-SCION-network-full.tsv"
    File network_rds = "out/network.rds"
    File target_rds = "out/target.rds"
    File reg_rds = "out/reg.rds"
    File cluster_assignment_rds = "out/cluster_assignment.rds"
    File params_rds = "out/params.rds"
  }

  runtime {
    container : "broadcptacdev/panoply_scion:latest"
    memory : "~{memory}GB"
    disks : "local-disk ~{disk_space} SSD"
    cpu : num_cores
    preemptible : num_preemptions
  }
}

task panoply_scion_run_permutation {
  input {
    File target_rds
    File reg_rds
    File cluster_assignment_rds
    File params_rds
    Int index
    Int base_seed
    String permute_dim
    Int num_preemptions
  }

  command <<<
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/scion_run_permutation.R \
      --target_rds ~{target_rds} \
      --reg_rds ~{reg_rds} \
      --cluster_assignment_rds ~{cluster_assignment_rds} \
      --params_rds ~{params_rds} \
      --index ~{index} \
      --base_seed ~{base_seed} \
      --permute_dim ~{permute_dim} \
      --num_cores 1 \
      --out_dir out
  >>>

  output {
    File permutation_rds = glob("out/permutation_*.rds")[0]
  }

  runtime {
    container : "broadcptacdev/panoply_scion:latest"
    memory : "8GB"
    disks : "local-disk 20 SSD"
    cpu : 1
    preemptible : num_preemptions
  }
}

task panoply_scion_aggregate_fdr {
  input {
    File network_rds
    Array[File] permutation_files
    Float target_fdr
    Int num_preemptions
  }

  command <<<
    set -euo pipefail
    mkdir -p permutations
    for f in ~{sep(" ", permutation_files)}; do
      cp "$f" permutations/
    done

    Rscript /prot/proteomics/Projects/PGDAC/src/scion_aggregate_fdr.R \
      --network_rds ~{network_rds} \
      --permutation_dir permutations \
      --target_fdr ~{target_fdr} \
      --out_dir out
  >>>

  output {
    File thresholded_network_tsv = "out/thresholded_network.tsv"
    File fdr_curve_png = "out/fdr_curve.png"
    File weight_comparison_png = "out/weight_comparison.png"
    File? network_plot_png = "out/network_plot.png"
    File fdr_result_rds = "out/fdr_result.rds"
  }

  runtime {
    container : "broadcptacdev/panoply_scion:latest"
    memory : "8GB"
    disks : "local-disk 20 SSD"
    cpu : 1
    preemptible : num_preemptions
  }
}
