#!/usr/bin/env Rscript
# Runs the real (unpermuted) SCION network and saves everything a downstream
# permutation shard needs: the processed target/regulator matrices, the fixed
# cluster assignment, and the exact params used (so permutations stay
# consistent with the real network -- see compute_fdr_threshold()).
#
# target/reg matrices are always read as GCT here and re-derived into plain
# TSVs via scion_prepare_gct.R first (see that file for why: PTM omes need
# their row names rebuilt from rdesc metadata, not used as raw GCT row IDs),
# then fed to SCION as format = "delimited".

suppressMessages({
  library(optparse)
  library(SCION)
})
# fixed deployment path (see panoply_scion.wdl) -- matches how this script
# itself is always invoked, rather than a fragile self-locating dirname trick
source("/prot/proteomics/Projects/PGDAC/src/scion_prepare_gct.R")

option_list <- list(
  make_option("--target_data_file", type = "character"),
  make_option("--reg_data_file", type = "character"),
  make_option("--ome", type = "character"),
  make_option("--gene_id_col", type = "character", default = "geneSymbol"),
  make_option("--ptm_type", type = "character", default = "SM"),
  make_option("--target_genes_file", type = "character", default = NULL),
  make_option("--reg_genes_file", type = "character", default = NULL),
  make_option("--gene_list_header", type = "logical", default = TRUE),
  make_option("--clustering_data_file", type = "character", default = NULL),
  make_option("--clustering_method", type = "character", default = "none"),
  make_option("--clustering_threshold", type = "double", default = 0.5),
  make_option("--clusters_file", type = "character", default = NULL),
  make_option("--connect_hubs", type = "logical", default = TRUE),
  make_option("--weightthreshold", type = "double", default = 0),
  make_option("--normalize", type = "logical", default = TRUE),
  make_option("--num_cores", type = "integer", default = 1),
  make_option("--ptm_sep", type = "character", default = "."),
  make_option("--seed", type = "integer", default = 2020),
  make_option("--nb_trees", type = "integer", default = 10000),
  make_option("--out_dir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# phosphoproteome/acetylome/ubiquitylome regulator matrices need their PTM
# site reconstructed from rdesc metadata; a plain proteome regulator matrix
# still gets its row names re-derived from gene_id_col (not assumed to
# already equal the raw GCT row ID), just without any site-splitting.
is_ptm <- opt$ome %in% c("phosphoproteome", "acetylome", "ubiquitylome")

target_mat <- prepare_scion_matrix(opt$target_data_file, role = "target")
reg_mat <- prepare_scion_matrix(opt$reg_data_file, role = "reg", gene_id_col = opt$gene_id_col,
                                 is_ptm = is_ptm, ptm_type = opt$ptm_type, ptm_sep = opt$ptm_sep)

# For PTM omes, reg_mat's row names are "<gene><ptm_sep><site>" composite
# strings, but PANOPLY's TF_file lists plain gene symbols (one per omics-wide
# TF, not one per site) -- it has no way to know sites in advance. SCION's
# own reg_genes_file filtering matches against the FULL row name, so a
# plain-gene-symbol list would never match a PTM row name. Filter here
# instead, against just the gene portion, and don't pass reg_genes_file on
# to run_scion() for the PTM case.
reg_genes_file_for_run_scion <- opt$reg_genes_file
if (is_ptm && !is.null(opt$reg_genes_file)) {
  reg_genes <- if (isTRUE(opt$gene_list_header)) {
    utils::read.delim(opt$reg_genes_file, header = TRUE)[[1]]
  } else {
    readLines(opt$reg_genes_file)
  }
  reg_gene_part <- vapply(strsplit(rownames(reg_mat), opt$ptm_sep, fixed = TRUE), `[`, character(1), 1)
  reg_mat <- reg_mat[reg_gene_part %in% reg_genes, , drop = FALSE]
  reg_genes_file_for_run_scion <- NULL
}

target_tsv <- file.path(opt$out_dir, "target_prepared.tsv")
reg_tsv <- file.path(opt$out_dir, "reg_prepared.tsv")
utils::write.table(target_mat, target_tsv, sep = "\t", col.names = NA, quote = FALSE)
utils::write.table(reg_mat, reg_tsv, sep = "\t", col.names = NA, quote = FALSE)

result <- run_scion(
  target_data_file = target_tsv,
  reg_data_file = reg_tsv,
  target_genes_file = opt$target_genes_file,
  reg_genes_file = reg_genes_file_for_run_scion,
  gene_list_header = opt$gene_list_header,
  clustering_data_file = opt$clustering_data_file,
  format = "delimited",
  clustering_method = opt$clustering_method,
  clustering_threshold = opt$clustering_threshold,
  clusters_file = opt$clusters_file,
  connect_hubs = opt$connect_hubs,
  weightthreshold = opt$weightthreshold,
  normalize = opt$normalize,
  num.cores = opt$num_cores,
  ptm_sep = opt$ptm_sep,
  seed = opt$seed,
  permute = FALSE,
  nb.trees = opt$nb_trees
)

write_scion_network(result$network, file.path(opt$out_dir, "network.tsv"))
saveRDS(result$network, file.path(opt$out_dir, "network.rds"))
saveRDS(result$target, file.path(opt$out_dir, "target.rds"))
saveRDS(result$reg, file.path(opt$out_dir, "reg.rds"))
saveRDS(result$cluster_assignment, file.path(opt$out_dir, "cluster_assignment.rds"))
saveRDS(c(result$params, list(nb.trees = opt$nb_trees)), file.path(opt$out_dir, "params.rds"))

cat("Done. Real network:", nrow(result$network), "edges.\n")
