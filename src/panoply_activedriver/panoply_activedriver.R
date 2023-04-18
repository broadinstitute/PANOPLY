library(pacman)

p_load(data.table)
p_load(stringr)
p_load(tidyr)
p_load(dplyr)
p_load(glue)
p_load(yaml)

p_load(cmapR)
p_load(seqinr)
p_load(maftools)

p_load(ActiveDriver)


args <- commandArgs(TRUE)

# args <- c(
#     "/Users/kpham/Desktop/mutation_ptm/luad-v1.3/luad-comb-v1.3-phosphoproteome-filtered-combined-bridging-batch-correct.gct",
#     "/Users/kpham/Desktop/mutation_ptm/luad-v1.3/luad-comb-v1.3-wxs_snv_maftools.maf",
#     "/Users/kpham/Desktop/mutation_ptm/Gencode.v34.pc_translations_clean3nr.fasta",
#     "/Users/kpham/Desktop/PANOPLae/PANOPLY_activedriver_mimp/src/panoply_common/master-parameters.yaml",
#     8
# )

# source("/Users/kpham/Desktop/PANOPLae/PANOPLY_activedriver_mimp/src/panoply_activedriver/activedriver_formatting.R")
source("/prot/proteomics/Projects/PGDAC/src/activedriver_formatting.R")


ptm_gct = as.character(args[1])
maf_file = as.character(args[2])
fasta_path = as.character(args[3])
yaml_file = as.character(args[4])
num_cores = as.integer(args[5])

yaml_params = read_yaml(yaml_file)

search_engine = yaml_params$panoply_activedriver$search_engine  # not currently used but might be helpful in the future
group_by_column = yaml_params$panoply_activedriver$group_by_column
residues = yaml_params$panoply_activedriver$residues
level = yaml_params$panoply_activedriver$level
mutation_AA_change_colname = yaml_params$panoply_activedriver$mutation_AA_change_colname
sample_id_col = yaml_params$panoply_activedriver$sample_id_col 
transcript_id_col = yaml_params$panoply_activedriver$transcript_id_col
other_transcript_col = yaml_params$panoply_activedriver$other_transcript_col
region_pval_adj_method = yaml_params$panoply_activedriver$region_pval_adj_method


# parse out all PTM sites for SpectrumMill or FragPipe generated GCTs for CPTAC
# assumes row IDs in format: "ENSP00000367263.4_S93s_1_1_93_93", "ENSP00000353114.4_T296tT298t_2_1_294_298"
ptm_sites <- get_ptm_sites(ptm_gct, id_column = group_by_column, residues = residues)

fasta <- read.fasta(fasta_path, seqtype = "AA", as.string = T, strip.desc = T, set.attributes = F)
transcript_map <- get_transcript_map(fasta)
sequences <- get_sequences(fasta, level = level)

mutations_all <- get_mutations_all(
    maf_file = maf_file,
    sample_col = sample_id_col,
    transcript_col = transcript_id_col,
    other_transcript_col = other_transcript_col
)

# ActiveDriver requires this column to be named "gene", even if we run analysis on protein-level
if (level != "transcript") {  # if not transcript then map to other levels ("protein" or "gene")
    mutations_all$gene <- transcript_map[mutations_all$gene, level]
}
mutations_all <- mutations_all %>% filter(!is.na(gene))

# TODO: allow input of disordered sequences
disordered_sequences <- get_disordered_regions_zero(sequences)

dir.create("activedriver_results_dir")
setwd("activedriver_results_dir")

results_all <- ActiveDriver(sequences, disordered_sequences, mutations_all, ptm_sites, mc.cores = num_cores)
results_all$merged_report$active_region_p_adj <- p.adjust(results_all$merged_report$active_region_p, method = region_pval_adj_method)
saveRDS(results_all, "results_all.rds")

protein_map <- get_protein_map(transcript_map)
ranked_genes <- AD_extract_top_genes(results_all, protein_map)
write.csv(ranked_genes, "top_genes_mutation_near_ptm.csv")