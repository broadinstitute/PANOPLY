# test_parameter_manager_R.r
### Commandline inputs:
library(optparse)

option_list <- list(
  # Every use:
  make_option(c("--module"), type = "character", dest = 'module', help = "The module or workflow name."),
  make_option(c("--master_yaml"), type="character", dest = 'master_yaml', help="master-parameters.yaml file to provide defaults"),
  # Global:
  make_option(c("--ndigits"), type = "integer", dest = 'ndigits', help = "Output precision for gct tables"),
  make_option(c("--na_max"), type = "double", dest = 'na_max', help = "maximum allowed NA values (per protein/site/row), can be fraction"),
  make_option(c("--sample_na_max"), type = "double", dest = 'sample_na_max', help = "maximum allowed fraction of NA values per sample/column; pipeline error if violated"),
  make_option(c("--min_numratio_fraction"), type = "double", dest = 'min_numratio_fraction', help = "fraction of samples in which min. numratio should be present to retain protein/phosphosite"),
  make_option(c("--nmiss_factor"), type = "double", dest = 'nmiss_factor', help = "for some situations, a more stringent condition is needed"),
  make_option(c("--sd_filter_threshold"), type = "double", dest = 'sd_filter_threshold', help = "SD threshold for SD filtering"),
  make_option(c("--duplicate_gene_policy"), type="character", dest = 'duplicate_gene_policy', help="Policy for combining/collapsing duplicate gene names"),
  make_option(c("--gene_id_col"), type="character", dest = 'gene_id_col', help="column name containing gene id"),
  make_option(c("--organism"), type="character", dest = 'organism', help="organism of study"),
  # All Command line possibles:
    # parse_SM_table:
  make_option(c("--label_type"), type="character", dest = 'label_type', help="label type for MS experiment. ex: TMT10", metavar="character"),
  make_option(c("--apply_sm_filter"), type="logical", dest = 'apply_sm_filter', help="if TRUE, apply numRatio based filter (use TRUE if input is SM ssv)"),
  make_option(c("--species_filter"), type="logical", dest = 'species_filter', help="if TRUE will discard all proteins/peptides not from homo sapiens"),
    #normalize_ms_data:
  make_option(c("--norm_method"), type = "character", dest = 'norm_method', help = "normalization method ex: '2comp', median, mean"),
  make_option(c("--alt_method"), type = "character", dest = 'alt_method', help = "alt.method for comparison -- filtered datasets not generated"),
    #rna_protein_correlation:
  make_option(c("--rna_sd_threshold"), type = "integer", dest = 'rna_sd_threshold', help = "for variation filter (set to NA to disable)"),
  make_option(c("--profile_plot_top_n"), type = "integer", dest = 'profile_plot_top_n', help = "defined in rna-seq-correlation.r"),
    #harmonize:
  make_option(c("--pome_gene_id_col"), type = "character", dest = 'pome_gene_id_col', help = "gene id column for harmonize pome data"),
  make_option(c("--cna_gene_id_col"), type = "character", dest = 'cna_gene_id_col', help = "gene id column for harmonize cna data"),
  make_option(c("--rna_gene_id_col"), type = "character", dest = 'rna_gene_id_col', help = "gene id column for harmonize rna data"),
    #association:
  make_option(c("--fdr_assoc"), type = "double", dest = 'fdr_assoc', help = "association module fdr value"),
    #association_report:
  make_option(c("--fdr_value"), type = "double", dest = 'fdr_value', help = "association report module fdr value"),
    #sample_qc:
  make_option(c("--cor_threshold"), type = "double", dest = 'cor_threshold', help = "correlation threshold for sample_qc"),
    #cna_analysis:
  make_option(c("--pe_max_default"), type = "integer", dest = 'pe_max_default', help = "default max precessors/jobs"),
  make_option(c("--fdr_cna_corr"), type = "double", dest = 'fdr_cna_corr', help = "fdr significance cutoff value set for cna analysis"),
  make_option(c("--min_cna_N"), type = "integer", dest = 'min_cna_N', help = "for cna-analysis-setup.r"),
    #cna_correlation_report:
  make_option(c("--cna_report_fdr"), type = "double", dest = 'cna_report_fdr', help = "fdr cutoff for cna reports."),
    #rna_protein_correlation_report:
  make_option(c("--rna_report_fdr"), type = "double", dest = 'rna_report_fdr', help = "fdr cutoff for rna reports."),
    #cons_cluster:
  make_option(c("--clustering_sd_threshold"), type = "integer", dest = "clustering_sd_threshold", help = "threshold for filtering data before consensus clustering"),
  make_option(c("--clustering_na_threshold"), type = "double", dest = "clustering_na_threshold", help = "max fraction of missing values for clustering; rest are imputed"),
  #immune_analysis:
  make_option(c("--immune_enrichment_fdr"), type = "double", dest = "immune_enrichment_fdr", help = "fdr value for immune analysis"),
  make_option(c("--immune_enrichment_subgroups"), type = "character", dest = "immune_enrichment_subgroups", help = "immune_enrichment_groups for immune analysis"),
  make_option(c("--immune_heatmap_width"), type = "integer", dest = "immune_heatmap_width", help = "immune_heatmap_width for immune analysis"),
  make_option(c("--immune_heatmap_height"), type = "integer", dest = "immune_heatmap_height", help = "immune_heatmap_height for immune analysis"),  
    #cmap_analysis:
  make_option(c("--cna_threshold"), type = "double", dest = 'cna_threshold', help = "copy number up/down if abs (log2(copy number) - 1) is > cna.threshold (default: 0.3)"),
  make_option(c("--cna_effects_threshold"), type = "integer", dest = 'cna_effects_threshold', help = "min number of samples with up/down copy number to include gene for CMAP analysis (default: 15)"),
  make_option(c("--min_sigevents"), type = "integer", dest = 'min_sigevents', help = "gene must have at least this many significant trans events to be considered (default: 20)"),
  make_option(c("--max_sigevents"), type = "integer", dest = 'max_sigevents', help = "if a gene has more then max.sigevents, the top max.sigevents will be chosen"),
  make_option(c("--top_N"), type = "integer", dest = 'top_N', help = "maximum number of genes to run CMAP enrichment on (default: 500)"),
  make_option(c("--fdr_pvalue"), type = "double", dest = 'fdr_pvalue', help = "FDR for CNA correlations (default: 0.05)"),
  make_option(c("--log_transform"), type = "logical", dest = 'log_transform', help = "if TRUE, log transform input data"),
  make_option(c("--must_include_genes"), dest = 'must_include_genes', help = "genes that must be included in the CMAP analysis (vector, default: NULL)"),
  make_option(c("--cis_fdr"), type = "character", dest = 'cis_fdr', help = "FDR for cis-correlations (default: 0.05)"),
  make_option(c("--legacy_score"), type = "logical", dest = 'legacy_score', help = "if TRUE, legacy connectivity score will be calculated (using mean rank points), with permutation FDR. if FALSE, enrichement will be based on fisher test of outlier scores, with BH-FDR (default)"),
  make_option(c("--rankpt_n"), type = "integer", dest = 'rankpt_n', help = "number of CMAP profiles to consider when calculating mean rank point (default: 4)"),
  make_option(c("--mean_rankpt_threshold"), type = "integer", dest = 'mean_rankpt_threshold', help = "min value of mean rank point for gene signature to be considered enriched (default: 85)"),
  make_option(c("--cmap_fdr"), type = "double", dest = 'cmap_fdr', help = "BH-FDR threshold for fisher test of outlier scores, for gene to be considered enriched"),
  make_option(c("--alpha"), type = "character", dest = 'alpha', help = "p-value threshold for cmap profile zscores and enrichments"),
  #mo_nmf:
  make_option(c("--mo_nmf_kmin"), type = "integer", dest = 'mo_nmf_kmin', help = "mo nmf clustering kmin value."),
  make_option(c("--mo_nmf_kmax"), type = "integer", dest = 'mo_nmf_kmax', help = "mo nmf clustering kmax value."),
  make_option(c("--mo_nmf_exclude_2"), type = "logical", dest = 'mo_nmf_exclude_2', help = "mo nmf clustering exclude_2 value."),
  make_option(c("--mo_nmf_core_membership"), type = "double", dest = 'mo_nmf_core_membership', help = "mo nmf clustering core_membership value."),
  make_option(c("--mo_nmf_nrun"), type = "integer", dest = 'mo_nmf_nrun', help = "mo nmf nrun value."),
  make_option(c("--mo_nmf_seed"), type = "character", dest = 'mo_nmf_seed', help = "mo nmf seed value."),
  make_option(c("--mo_nmf_method"), type = "character", dest = 'mo_nmf_method', help = "mo nmf method value."),
  make_option(c("--mo_nmf_bnmf"), type = "logical", dest = 'mo_nmf_bnmf', help = "mo nmf bmnf value."),
  make_option(c("--mo_nmf_cat_anno"), dest = 'mo_nmf_cat_anno', help = "mo nmf annotations to include in heatmaps and enrichments, categorical."),
  make_option(c("--mo_nmf_cont_anno"), dest = 'mo_nmf_cont_anno', help = "mo nmf annotations to include in heatmaps and enrichments, continuous."),
  make_option(c("--mo_nmf_cat_colors_cat1_level1"), type = "character", dest = 'mo_nmf_cat_colors_cat1_level1', help = "mo nmf colors for cat_anno."),
  make_option(c("--mo_nmf_cat_colors_cat1_level2"), type = "character", dest = 'mo_nmf_cat_colors_cat1_level2', help = "mo nmf colors for cat_anno."),
  make_option(c("--mo_nmf_cat_colors_cat2_level1"), type = "character", dest = 'mo_nmf_cat_colors_cat2_level1', help = "mo nmf colors for cat_anno."),
  make_option(c("--mo_nmf_cat_colors_cat2_level2"), type = "character", dest = 'mo_nmf_cat_colors_cat2_level2', help = "mo nmf colors for cat_anno."),
  make_option(c("--mo_nmf_hm_cw"), type = "integer", dest = 'mo_nmf_hm_cw', help = "mo nmf heatmap cell width value."),
  make_option(c("--mo_nmf_hm_ch"), type = "integer", dest = 'mo_nmf_hm_ch', help = "mo nmf heatmap cell height value."),
  make_option(c("--mo_nmf_hm_max_val"), type = "integer", dest = 'mo_nmf_hm_max_val', help = "mo nmf heatmap max value at which the data will be capped in heatmap."),
  make_option(c("--mo_nmf_hm_max_val_z"), type = "integer", dest = 'mo_nmf_hm_max_val_z', help = "mo nmf heatmap max value if z-scored."),
  make_option(c("--mo_nmf_filt_mode"), type = "character", dest = 'mo_nmf_filt_mode', help = "mo nmf filtering mode value."),
  make_option(c("--mo_nmf_sd_filt"), type = "double", dest = 'mo_nmf_sd_filt', help = "mo nmf filtering mode value."),
  make_option(c("--mo_nmf_z_score"), type = "logical", dest = 'mo_nmf_z_score', help = "mo nmf filtering z score value."),
  make_option(c("--mo_nmf_impute"), type = "logical", dest = 'mo_nmf_impute', help = "mo nmf imputation value."),
  make_option(c("--mo_nmf_impute_k"), type = "integer", dest = 'mo_nmf_impute_k', help = "mo nmf imputation k value."),
  make_option(c("--mo_nmf_max_na_row"), type = "double", dest = 'mo_nmf_max_na_row', help = "mo nmf imputation max na row value."),
  make_option(c("--mo_nmf_max_na_col"), type = "double", dest = 'mo_nmf_max_na_col', help = "mo nmf imputation max na col value."),
  make_option(c("--mo_nmf_nmf_only"), type = "logical", dest = 'mo_nmf_nmf_only', help = "mo nmf value if TRUE, no downstream analysis of clustering results will be performed."),
  # ssGSEA_projection:
  make_option(c("--ssgsea_nperm"), type = "integer", dest = 'ssgsea_nperm', help = "ssgsea projection nperm value."),
  make_option(c("--ssgsea_weight"), type = "double", dest = 'ssgsea_weight', help = "ssgsea projection weight value."),
  make_option(c("--ssgsea_norm"), type = "character", dest = 'ssgsea_norm', help = "ssgsea projection normalization type (if any)."),
  make_option(c("--ssgsea_correl_type"), type = "character", dest = 'ssgsea_correl_type', help = "ssgsea projection normalization type (if any)."),
  make_option(c("--ssgsea_statistic"), type = "character", dest = 'ssgsea_statistic', help = "ssgsea statistic method (ie: area.under.RES, Kolmogorov-Smirnov ."),
  make_option(c("--ssgsea_output_score_types"), type = "character", dest = 'ssgsea_output_score_types', help = "ssgsea score type."),
  make_option(c("--ssgsea_min_overlap"), type = "integer", dest = 'ssgsea_min_overlap', help = "ssgsea Minimal overlap between signature and data set."),
  make_option(c("--ssgsea_extended_output"), type = "logical", dest = 'ssgsea_extended_output', help = "ssgsea If TRUE additional stats on signature coverage etc. will be included as row annotations in the GCT results files."),
  make_option(c("--ssgsea_export_signat_gct"), type = "logical", dest = 'ssgsea_export_signat_gct', help = "ssgsea For each signature export expression GCT files."),
  make_option(c("--ssgsea_global_fdr"), type = "logical", dest = 'ssgsea_global_fdr', help = "ssgsea If TRUE global FDR across all data columns is calculated."),
  make_option(c("--ssgsea_multi_core"), type = "logical", dest = 'ssgsea_multi_core', help = "ssgsea If TRUE processing will be parallelized across gene sets using (N-1) CPU cores"),
  # blacksheep_module:
  make_option(c("--blacksheep_apply_filtering"), type = "logical", dest = 'blacksheep_apply_filtering', help = "blacksheep_apply_filtering"),
  make_option(c("--blacksheep_identifiers_file"), dest = 'blacksheep_identifiers_file', help = "blacksheep_identifiers_file"),
  make_option(c("--blacksheep_groups_file"), dest = 'blacksheep_groups_file', help = "blacksheep_groups_file"),
  make_option(c("--blacksheep_fraction_samples_cutoff"), type = "double", dest = 'blacksheep_fraction_samples_cutoff', help = "blacksheep_fraction_samples_cutoff"),
  make_option(c("--blacksheep_fdr_value"), type = "double", dest = 'blacksheep_fdr_value', help = "blacksheep_fdr_value")
  )


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

### Package to read/write YAML parameter file:
if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load('yaml')

# MODULE NAMES:
# pipeline --evals all parameters at once
# parse_sm_table
# normalize_ms_data
# rna_protein_correlation
# harmonize
# sample_qc
# cna_analysis
# cna_corr_report
# rna_corr_report
# cons_clust
# association
# association_report
# cmap_analysis
# immune_analysis
# blacksheep
# mo_nmf
# ssgsea_projection

### FUNCTIONS:

### HELPER FUNCTIONS:
# all values are optional other than help, module, and master_yaml so once those are removed
# opt will be of length 0 if no command line overwrites were given!
check_if_any_command_line <- function(opt){
  opt_copy <- opt
  opt_copy$help <- NULL
  opt_copy$module <- NULL
  opt_copy$master_yaml <- NULL
  if (length(opt_copy) == 0){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

# Read in the yaml file:
read_yaml_file <- function(yaml){
  open_yaml <- read_yaml(yaml)
  return(open_yaml)
}

# Write out final yaml file:
write_yaml_file <- function(yaml){
  yaml_name <- file('final_output_params.yaml', "w")
  write_yaml(yaml, yaml_name, handlers = list(
    logical = function(x) {
      result <- ifelse(x, "TRUE", "FALSE")
      class(result) <- "verbatim"
      return(result)
    }
  ))
  close(yaml_name)
}

######################################## CHECK PARAMS FUNCTIONS ################################
#Checks if any of the global vals are being changed from commandline. if yes will change them in yaml if no returns org yaml
check_global_params <- function(opt,yaml){
  if (!is.null(opt$ndigits) | !is.null(opt$na_max) | !is.null(opt$sample_na_max) | !is.null(opt$min_numratio_fraction) | !is.null(opt$nmiss_factor) | !is.null(opt$sd_filter_threshold) | !is.null(opt$duplicate_gene_policy | !is.null(opt$gene_id_col) | !is.null(opt$organism))){
    yaml <- change_global_params(opt,yaml)
    return(yaml)
  }else{
    return(yaml)
  }
}

#Is called when globals are being changed via command line. Will change any globals and return updated yaml
change_global_params <- function(opt,yaml){
  if (!is.null(opt$ndigits)){
    yaml$global_parameters$output_precision$ndigits <- opt$ndigits
  }
  if (!is.null(opt$na_max)){
    yaml$global_parameters$missing_values_and_filtering$na_max <- opt$na_max
  }
  if (!is.null(opt$sample_na_max)){
    yaml$global_parameters$missing_values_and_filtering$sample_na_max <- opt$sample_na_max
  }
  if (!is.null(opt$min_numratio_fraction)){
    yaml$global_parameters$missing_values_and_filtering$min_numratio_fraction <- opt$min_numratio_fraction
  }
  if (!is.null(opt$nmiss_factor)){
    yaml$global_parameters$missing_values_and_filtering$nmiss_factor <- opt$nmiss_factor
  }
  if (!is.null(opt$sd_filter_threshold)){
    yaml$global_parameters$missing_values_and_filtering$sd_filter_threshold <- opt$sd_filter_threshold
  }
  if (!is.null(opt$duplicate_gene_policy)){
    yaml$global_parameters$gene_mapping$duplicate_gene_policy <- opt$duplicate_gene_policy
  }
  if (!is.null(opt$gene_id_col)){
    yaml$global_parameters$gene_mapping$gene_id_col <- opt$gene_id_col
  }
  if (!is.null(opt$organism)){
    yaml$global_parameters$organism <- opt$organism
  }
  return(yaml)
}

# parse_sm_table:
check_parse_sm_params <- function(opt, yaml){
  if (!is.null(opt$label_type) | !is.null(opt$apply_sm_filter) | !is.null(opt$species_filter)){
    if (!is.null(opt$label_type)){
      yaml$panoply_parse_sm_table$label_type_for_ms_exp$label_type <- opt$label_type
    }
    if (!is.null(opt$apply_sm_filter)) {
      yaml$panoply_parse_sm_table$apply_sm_filter <- opt$apply_sm_filter
    }
    if (!is.null(opt$species_filter)) {
      yaml$panoply_parse_sm_table$species_filter <- opt$species_filter
    }
    return(yaml)
  }else{
    return(yaml)
  }
}

# normalize_sm_data:
check_normalize_sm_params <- function(opt, yaml){
  if (!is.null(opt$norm_method) | !is.null(opt$alt_method)){
    if (!is.null(opt$norm_method)){
      yaml$panoply_normalize_ms_data$normalization$norm_method <- opt$norm_method
    }
    if (!is.null(opt$alt_method)) {
      yaml$panoply_normalize_ms_data$normalization$alt_method <- opt$alt_method
    }
    return(yaml)
  }else{
    return(yaml)
  }
}

# rna_protein_correlation:
check_rna_protein_correlation_params <- function(opt, yaml){
  if (!is.null(opt$rna_sd_threshold) | !is.null(opt$profile_plot_top_n)){
    if (!is.null(opt$rna_sd_threshold)){
      yaml$panoply_rna_protein_correlation$rna$rna_sd_threshold <- opt$rna_sd_threshold
    }
    if (!is.null(opt$profile_plot_top_n)) {
      yaml$panoply_rna_protein_correlation$rna$profile_plot_top_n <- opt$profile_plot_top_n
    }
    return(yaml)
  }else{
    return(yaml)
  }
}

# harmonize:
check_harmonize_params <- function(opt, yaml){
  if (!is.null(opt$pome_gene_id_col) | !is.null(opt$cna_gene_id_col) | !is.null(opt$rna_gene_id_col)){
    if (!is.null(opt$pome_gene_id_col)){
      yaml$panoply_harmonize$pome_gene_id_col <- opt$pome_gene_id_col
    }
    if (!is.null(opt$cna_gene_id_col)) {
      yaml$panoply_harmonize$cna_gene_id_col <- opt$cna_gene_id_col
    }
    if (!is.null(opt$rna_gene_id_col)) {
      yaml$panoply_harmonize$rna_gene_id_col <- opt$rna_gene_id_col
    }
    return(yaml)
  }else{
    return(yaml)
  }
}

# sample_qc:
check_sample_qc_params <- function(opt, yaml){
  if (!is.null(opt$cor_threshold)){
    yaml$panoply_sample_qc$cor_threshold <- opt$cor_threshold
    return(yaml)
  }else{
    return(yaml)
  }
}

# cna_analysis
check_cna_analysis_params <- function(opt, yaml){
  if (!is.null(opt$pe_max_default) | !is.null(opt$fdr_cna_corr) | !is.null(opt$min_cna_N)){
    if (!is.null(opt$pe_max_default)){
      yaml$panoply_cna_analysis$cna_parallelism$pe_max_default <- opt$pe_max_default
    }
    if (!is.null(opt$fdr_cna_corr)) {
      yaml$panoply_cna_analysis$fdr_cna_corr <- opt$fdr_cna_corr
    }
    if (!is.null(opt$min_cna_N)) {
      yaml$panoply_cna_analysis$min_cna_N <- opt$min_cna_N
    }
    return(yaml)
  }else{
    return(yaml)
  }
}

# cna_corr_report
check_cna_corr_report_params <- function(opt, yaml){
  if (!is.null(opt$cna_report_fdr)){
    yaml$panoply_cna_correlation_report$fdr <- opt$cna_report_fdr
    return(yaml)
  }else{
    return(yaml)
  }
}

# rna_corr_report
check_rna_corr_report_params <- function(opt, yaml){
  if (!is.null(opt$rna_report_fdr)){
    yaml$panoply_rna_protein_correlation_report$fdr <- opt$rna_report_fdr
    return(yaml)
  }else{
    return(yaml)
  }
}

# cons_clust:
check_cons_clust_params <- function(opt, yaml){
  if (!is.null(opt$clustering_sd_threshold) | !is.null(opt$clustering_na_threshold)){
    if (!is.null(opt$clustering_sd_threshold)){
      yaml$panoply_cons_cluster$clustering_sd_threshold <- opt$clustering_sd_threshold
    }
    if (!is.null(opt$clustering_na_threshold)) {
      yaml$panoply_cons_cluster$clustering_na_threshold <- opt$clustering_na_threshold
    }
    return(yaml)
  }else{
    return(yaml)
  }
}

# cmap_analysis:
check_cmap_analysis_params <- function(opt, yaml){
  if (!is.null(opt$cna_threshold) | !is.null(opt$cna_effects_threshold) | !is.null(opt$min_sigevents) | !is.null(opt$max_sigevents) | !is.null(opt$top_N) | !is.null(opt$fdr_pvalue) | !is.null(opt$log_transform) | !is.null(opt$must_include_genes) | !is.null(opt$cis_fdr) | !is.null(opt$legacy_score) | !is.null(opt$rankpt_n) | !is.null(opt$mean_rankpt_threshold) | !is.null(opt$cmap_fdr) | !is.null(opt$alpha)){
    if (!is.null(opt$cna_threshold)){
      yaml$panoply_cmap_analysis$cna_threshold <- opt$cna_threshold
    }
    if (!is.null(opt$cna_effects_threshold)){
      yaml$panoply_cmap_analysis$cna_effects_threshold <- opt$cna_effects_threshold
    }
    if (!is.null(opt$min_sigevents)){
      yaml$panoply_cmap_analysis$min_sigevents <- opt$min_sigevents
    }
    if (!is.null(opt$max_sigevents)){
      yaml$panoply_cmap_analysis$max_sigevents <- opt$max_sigevents
    }
    if (!is.null(opt$top_N)){
      yaml$panoply_cmap_analysis$top_N <- opt$top_N
    }
    if (!is.null(opt$fdr_pvalue)){
      yaml$panoply_cmap_analysis$fdr_pvalue <- opt$fdr_pvalue
    }
    if (!is.null(opt$log_transform)){
      yaml$panoply_cmap_analysis$log_transform <- opt$log_transform
    }
    if (!is.null(opt$must_include_genes)){
      yaml$panoply_cmap_analysis$must_include_genes <- opt$must_include_genes
    }
    if (!is.null(opt$cis_fdr)){
      yaml$panoply_cmap_analysis$cis_fdr <- opt$cis_fdr
    }
    if (!is.null(opt$legacy_score)){
      yaml$panoply_cmap_analysis$legacy_score <- opt$legacy_score
    }
    if (!is.null(opt$rankpt_n)){
      yaml$panoply_cmap_analysis$rankpt_n <- opt$rankpt_n
    }
    if (!is.null(opt$mean_rankpt_threshold)){
      yaml$panoply_cmap_analysis$mean_rankpt_threshold <- opt$mean_rankpt_threshold
    }
    if (!is.null(opt$cmap_fdr)){
      yaml$panoply_cmap_analysis$cmap_fdr <- opt$cmap_fdr
    }
    if (!is.null(opt$alpha)){
      yaml$panoply_cmap_analysis$alpha <- opt$alpha
    }
    return(yaml)
  }else{
    return(yaml)
  }
  
}

# immune_analysis:
check_immune_analysis_params <- function(opt, yaml){
  if (!is.null(opt$immune_enrichment_fdr) | !is.null(opt$immune_enrichment_subgroups) | !is.null(opt$immune_heatmap_width) | !is.null(opt$immune_heatmap_height)){
    if (!is.null(opt$immune_enrichment_fdr)){
      yaml$panoply_immune_analysis$immune_enrichment_fdr <- opt$immune_enrichment_fdr
    }
    if (!is.null(opt$immune_enrichment_subgroups)) {
      yaml$panoply_immune_analysis$immune_enrichment_subgroups <- opt$immune_enrichment_subgroups
    }
    if (!is.null(opt$immune_heatmap_width)) {
      yaml$panoply_immune_analysis$immune_heatmap_width <- opt$immune_heatmap_width
    }
    if (!is.null(opt$immune_heatmap_height)) {
      yaml$panoply_immune_analysis$immune_heatmap_height <- opt$immune_heatmap_height
    }
    return(yaml)
  }else{
    return(yaml)
  }
}

# mo_nmf:
check_mo_nmf_params <- function(opt, yaml){
  if (!is.null(opt$mo_nmf_kmin)){
    yaml$panoply_mo_nmf$kmin <- opt$mo_nmf_kmin
  }
  if (!is.null(opt$mo_nmf_kmax)){
    yaml$panoply_mo_nmf$kmax <- opt$mo_nmf_kmax
  }
  if (!is.null(opt$mo_nmf_exclude_2)){
    yaml$panoply_mo_nmf$exclude_2 <- opt$mo_nmf_exclude_2
  }
  if (!is.null(opt$mo_nmf_core_membership)){
    yaml$panoply_mo_nmf$core_membership <- opt$mo_nmf_core_membership
  }
  if (!is.null(opt$mo_nmf_nrun)){
    yaml$panoply_mo_nmf$nrun <- opt$mo_nmf_nrun
  }
  if (!is.null(opt$mo_nmf_seed)){
    yaml$panoply_mo_nmf$seed <- opt$mo_nmf_seed
  }
  if (!is.null(opt$mo_nmf_method)){
    yaml$panoply_mo_nmf$method <- opt$mo_nmf_method
  }
  if (!is.null(opt$mo_nmf_bnmf)){
    yaml$panoply_mo_nmf$bnmf <- opt$mo_nmf_bnmf
  }
  if (!is.null(opt$mo_nmf_cat_anno)){
    yaml$panoply_mo_nmf$cat_anno <- opt$mo_nmf_cat_anno
  }
  if (!is.null(opt$mo_nmf_cont_anno)){
    yaml$panoply_mo_nmf$cont_anno <- opt$mo_nmf_cont_anno
  }
  if (!is.null(opt$mo_nmf_cat_colors_cat1_level1)){
    yaml$panoply_mo_nmf$cat_colors$cat1$level1 <- opt$mo_nmf_cat_colors_cat1_level1
  }
  if (!is.null(opt$mo_nmf_cat_colors_cat1_level2)){
    yaml$panoply_mo_nmf$cat_colors$cat1$level2 <- opt$mo_nmf_cat_colors_cat1_level2
  }
  if (!is.null(opt$mo_nmf_cat_colors_cat2_level1)){
    yaml$panoply_mo_nmf$cat_colors$cat2$level1 <- opt$mo_nmf_cat_colors_cat2_level1
  }
  if (!is.null(opt$mo_nmf_cat_colors_cat2_level2)){
    yaml$panoply_mo_nmf$cat_colors$cat2$level2 <- opt$mo_nmf_cat_colors_cat2_level2
  }
  if (!is.null(opt$mo_nmf_hm_cw)){
    yaml$panoply_mo_nmf$hm_cw <- opt$mo_nmf_hm_cw
  }
  if (!is.null(opt$mo_nmf_hm_ch)){
    yaml$panoply_mo_nmf$hm_ch <- opt$mo_nmf_hm_ch
  }
  if (!is.null(opt$mo_nmf_hm_max_val)){
    yaml$panoply_mo_nmf$hm_max_val <- opt$mo_nmf_hm_max_val
  }
  if (!is.null(opt$mo_nmf_hm_max_val_z)){
    yaml$panoply_mo_nmf$hm_max_val_z <- opt$mo_nmf_hm_max_val_z
  }
  if (!is.null(opt$mo_nmf_filt_mode)){
    yaml$panoply_mo_nmf$filt_mode <- opt$mo_nmf_filt_mode
  }
  if (!is.null(opt$mo_nmf_sd_filt)){
    yaml$panoply_mo_nmf$sd_filt <- opt$mo_nmf_sd_filt
  }
  if (!is.null(opt$mo_nmf_z_score)){
    yaml$panoply_mo_nmf$z_score <- opt$mo_nmf_z_score
  }
  if (!is.null(opt$mo_nmf_impute)){
    yaml$panoply_mo_nmf$impute <- opt$mo_nmf_impute
  }
  if (!is.null(opt$mo_nmf_impute_k)){
    yaml$panoply_mo_nmf$impute_k <- opt$mo_nmf_impute_k
  }
  if (!is.null(opt$mo_nmf_max_na_row)){
    yaml$panoply_mo_nmf$max_na_row <- opt$mo_nmf_max_na_row
  }
  if (!is.null(opt$mo_nmf_max_na_col)){
    yaml$panoply_mo_nmf$max_na_col <- opt$mo_nmf_max_na_col
  }
  if (!is.null(opt$mo_nmf_nmf_only)){
    yaml$panoply_mo_nmf$nmf_only <- opt$mo_nmf_nmf_only
  }
  return(yaml)
}

# ssgsea_projection:
check_ssgsea_projection_params <- function(opt, yaml){
  if (!is.null(opt$ssgsea_nperm)){
    yaml$panoply_ssgsea_projection$nperm <- opt$ssgsea_nperm
  }
  if (!is.null(opt$ssgsea_weight)){
    yaml$panoply_ssgsea_projection$weight <- opt$ssgsea_weight
  }
  if (!is.null(opt$ssgsea_norm)){
    yaml$panoply_ssgsea_projection$norm <- opt$ssgsea_norm
  }
  if (!is.null(opt$ssgsea_correl_type)){
    yaml$panoply_ssgsea_projection$correl_type <- opt$ssgsea_correl_type
  }
  if (!is.null(opt$ssgsea_statistic)){
    yaml$panoply_ssgsea_projection$statistic <- opt$ssgsea_statistic
  }
  if (!is.null(opt$ssgsea_output_score_type)){
    yaml$panoply_ssgsea_projection$output_score_type <- opt$ssgsea_output_score_type
  }
  if (!is.null(opt$ssgsea_min_overlap)){
    yaml$panoply_ssgsea_projection$min_overlap <- opt$ssgsea_min_overlap
  }
  if (!is.null(opt$ssgsea_extended_output)){
    yaml$panoply_ssgsea_projection$extended_output <- opt$ssgsea_extended_output
  }
  if (!is.null(opt$ssgsea_export_signat_gct)){
    yaml$panoply_ssgsea_projection$export_signat_gct <- opt$ssgsea_export_signat_gct
  }
  if (!is.null(opt$ssgsea_global_fdr)){
    yaml$panoply_ssgsea_projection$global_fdr <- opt$ssgsea_global_fdr
  }
  if (!is.null(opt$ssgsea_multi_core)){
    yaml$panoply_ssgsea_projection$multi_core <- opt$ssgsea_multi_core
  }
  return(yaml)
}

# blacksheep:
check_blacksheep_params <- function(opt, yaml){
  if (!is.null(opt$blacksheep_apply_filtering)){
    yaml$panoply_blacksheep$apply_filtering <- opt$blacksheep_apply_filtering
  }
  if (!is.null(opt$blacksheep_identifiers_file)){
    yaml$panoply_blacksheep$identifiers_file <- opt$blacksheep_identifiers_file
  }
  if (!is.null(opt$blacksheep_groups_file)){
    yaml$panoply_blacksheep$groups_file <- opt$blacksheep_groups_file
  }
  if (!is.null(opt$blacksheep_fraction_samples_cutoff)){
    yaml$panoply_blacksheep$fraction_samples_cutoff <- opt$blacksheep_fraction_samples_cutoff
  }
  if (!is.null(opt$blacksheep_fdr_value)){
    yaml$panoply_blacksheep$fdr_value <- opt$blacksheep_fdr_value
  }
  return(yaml)
}

# association:
check_association_params <- function(opt, yaml){
  if (!is.null(opt$fdr_assoc)){
    yaml$panoply_association$fdr_assoc <- opt$fdr_assoc
  }
  return(yaml)
}

#association_report:
check_association_report_params <- function(opt, yaml){
  if (!is.null(opt$fdr_value)){
    yaml$panoply_association_report$fdr_value <- opt$fdr_value
  }
  return(yaml)
}

# Checks all parameters (maybe use for final output yaml for whole pipeline?)
check_pipeline_params <- function(opt,yaml){
  yaml <- check_global_params(opt, yaml)
  yaml <- check_parse_sm_params(opt,yaml)
  yaml <- check_normalize_sm_params(opt,yaml)
  yaml <- check_rna_protein_correlation_params(opt,yaml)
  yaml <- check_harmonize_params(opt,yaml)
  yaml <- check_association_params(opt,yaml)
  yaml <- check_association_report_params(opt,yaml)
  yaml <- check_sample_qc_params(opt,yaml)
  yaml <- check_cna_analysis_params(opt,yaml)
  yaml <- check_cna_corr_report_params(opt,yaml)
  yaml <- check_rna_corr_report_params(opt,yaml)
  yaml <- check_cons_clust_params(opt,yaml)
  yaml <- check_cmap_analysis_params(opt,yaml)
  yaml <- check_immune_analysis_params(opt, yaml)
  yaml <- check_mo_nmf_params(opt, yaml)
  yaml <- check_ssgsea_projection_params(opt, yaml)
  yaml <- check_blacksheep_params(opt, yaml)
  return(yaml)
}

##################################################################################################
# Write to config functions:
# For all modules other than cmap:
write_custom_config <- function(yaml){
  #custom_config_path <- '/prot/proteomics/Projects/PGDAC/src/'
  #custom_config_path <- '/output/'
  output <- paste(paste('ndigits', '<-', yaml$global_parameters$output_precision$ndigits),
                paste('na.max', '<-', yaml$global_parameters$missing_values_and_filtering$na_max), 
                paste('sample.na.max', '<-', yaml$global_parameters$missing_values_and_filtering$sample_na_max),
                paste('min.numratio.fraction', '<-', yaml$global_parameters$missing_values_and_filtering$min_numratio_fraction),
                paste('nmiss.factor', '<-', yaml$global_parameters$missing_values_and_filtering$nmiss_factor),
                paste('sd.filter.threshold', '<-', yaml$global_parameters$missing_values_and_filtering$sd_filter_threshold),
                paste('duplicate.gene.policy', '<-', paste('"', yaml$global_parameters$gene_mapping$duplicate_gene_policy, '"', sep = '')),
                paste('gene.id.col', '<-', paste('"', yaml$global_parameters$gene_mapping$gene_id_col, '"', sep = '')),
                paste('organism', '<-', paste('"', yaml$global_parameters$organism, '"', sep = '')),
                #parse_sm_table:
                paste('label.type', '<-', paste('"', yaml$panoply_parse_sm_table$label_type_for_ms_exp$label_type, '"', sep = '')),
                'set.label.type (label.type)', #This needs to be run after label.type in config.r!
                paste('apply.sm.filter', '<-', yaml$panoply_parse_sm_table$apply_sm_filter),
                paste('species.filter', '<-', yaml$panoply_parse_sm_table$species_filter),
                #normalize_sm_data:
                paste('norm.method', '<-', paste('"', yaml$panoply_normalize_ms_data$normalization$norm_method, '"', sep = '')),
                paste('alt.method', '<-', paste('"', yaml$panoply_normalize_ms_data$normalization$alt_method, '"', sep = '')),
                'if (norm.method == alt.method) alt.method <- NULL', #This should be run after norm and alt methods are defined
                #rna_protein_correlation:
                paste('rna.sd.threshold', '<-', yaml$panoply_rna_protein_correlation$rna$rna_sd_threshold),
                paste('profile.plot.top.n', '<-', yaml$panoply_rna_protein_correlation$rna$profile_plot_top_n),
                #harmonize:
                paste('pome.gene.id.col', '<-', paste('"', yaml$panoply_harmonize$pome_gene_id_col, '"', sep = '')),
                paste('cna.gene.id.col', '<-', paste('"', yaml$panoply_harmonize$cna_gene_id_col, '"', sep = '')),
                paste('rna.gene.id.col', '<-', paste('"', yaml$panoply_harmonize$rna_gene_id_col, '"', sep = '')),
                #association:
                paste('assoc.fdr', '<-', yaml$panoply_association$fdr_assoc),
                #sample_qc:
                paste('cor.threshold', '<-', yaml$panoply_sample_qc$cor_threshold),
                #cna_analysis:
                paste('pe.max.default', '<-', yaml$panoply_cna_analysis$cna_parallelism$pe_max_default),
                paste('fdr_cna_corr', '<-', yaml$panoply_cna_analysis$fdr_cna_corr),
                paste('min.cna.N', '<-', yaml$panoply_cna_analysis$min_cna_N),
                #cna_correlation_report:
                paste('cna.report.fdr', '<-', yaml$panoply_cna_correlation_report$fdr),
                #rna_protein_correlation_report:
                paste('rna.report.fdr', '<-', yaml$panoply_rna_protein_correlation_report$fdr),
                #cons_cluster:
                paste('clustering.sd.threshold', '<-', yaml$panoply_cons_cluster$clustering_sd_threshold),
                paste('clustering.na.threshold', '<-', yaml$panoply_cons_cluster$clustering_na_threshold),
                #cmap_analysis: HAS ITS OWN FUNCTION below
                #immune_analysis:
                paste('immune.enrichment.fdr', '<-', yaml$panoply_immune_analysis$immune_enrichment_fdr),
                paste('immune.enrichment.subgroups', '<-', ifelse(is.null(yaml$panoply_immune_analysis$immune_enrichment_subgroups),'NULL', paste('"',yaml$panoply_immune_analysis$immune_enrichment_subgroups,'"', sep = ''))),
                paste('immune.heatmap.width', '<-', yaml$panoply_immune_analysis$immune_heatmap_width),
                paste('immune.heatmap.height', '<-', yaml$panoply_immune_analysis$immune_heatmap_height),
                #DEV_sample_annotations:
                paste('sample.id.col.name', '<-', paste('"',yaml$DEV_sample_annotation$sample_id_col_name, '"', sep = '')),
                paste('experiment.col.name', '<-', paste('"',yaml$DEV_sample_annotation$experiment_col_name, '"', sep = '')),
                paste('channel.col.name', '<-', paste('"',yaml$DEV_sample_annotation$channel_col_name, '"', sep = '')),
                paste('participant.col.name', '<-', paste('"',yaml$DEV_sample_annotation$participant_col_name, '"', sep = '')),
                paste('type.col.name', '<-', paste('"',yaml$DEV_sample_annotation$type_col_name, '"', sep = '')),
                paste('sampleQC.cls', '<-', yaml$DEV_sample_annotation$QC$sampleQC_cls),
                paste('qc.col', '<-', paste('"', yaml$DEV_sample_annotation$QC$qc_col, '"', sep = '')),
                paste('qc.pass.label', '<-', paste('"', yaml$DEV_sample_annotation$QC$qc_pass_label, '"', sep = '')),
                "if (type == 'proteome') {",
                paste('  id.col', '<-', yaml$DEV_sample_annotation$gct_file_ids$proteome$id_col),
                paste('  desc.col', '<-', yaml$DEV_sample_annotation$gct_file_ids$proteome$desc_col),
                paste('  min.numratio', '<-', yaml$DEV_sample_annotation$gct_file_ids$proteome$min_numratio),
                "} else {",
                paste('  id.col', '<-', yaml$DEV_sample_annotation$gct_file_ids$phosphoproteome$id_col),
                paste('  desc.col', '<-', yaml$DEV_sample_annotation$gct_file_ids$phosphoproteome$desc_col),
                paste('  min.numratio', '<-', yaml$DEV_sample_annotation$gct_file_ids$phosphoproteome$min_numratio),
                '}',
                #DEV_directory_and_file_names:
                paste('rna.output.prefix', '<-', paste('"', yaml$DEV_directory_and_file_names$rna_output_prefix, '"', sep = '')),
                paste('official.genesyms', '<-', paste('"', yaml$DEV_directory_and_file_names$gene_mapping$official_genesyms, '"', sep = '')),
                paste('protein.gene.map', '<-', paste('"', yaml$DEV_directory_and_file_names$gene_mapping$protein_gene_map, '"', sep = '')),
                paste('data.dir', '<-', paste("'",'../', yaml$DEV_directory_and_file_names$data_dir,"'", sep = '')),
                paste('parse.dir', '<-', paste("'",'../', yaml$DEV_directory_and_file_names$parse_dir,"'", sep = '')),
                paste('norm.dir', '<-', paste("'",'../', yaml$DEV_directory_and_file_names$norm_dir,"'", sep = '')),
                paste('harmonize.dir', '<-', paste("'",'../', yaml$DEV_directory_and_file_names$harmonize_dir,"'", sep = '')),
                paste('expt.design.file <- file.path (data.dir, ', paste("'", yaml$DEV_directory_and_file_names$expt_design_file, "'", sep = ''), ')', sep = ''),
                paste('rna.data.file <- file.path (data.dir, ', paste("'", yaml$DEV_directory_and_file_names$rna_data_file, "'", sep = ''), ')', sep = ''),
                paste('cna.data.file <- file.path (data.dir, ', paste("'", yaml$DEV_directory_and_file_names$cna_data_file, "'", sep = ''), ')', sep = ''),
                sep = "\n")
  write(output, 'config-custom.r')
}

# Write parameters to cmap config:
write_custom_cmap_config <- function(yaml){
  #custom_config_path <- '/prot/proteomics/Projects/PGDAC/src/'
  #custom_config_path <- '~/output/'
  output <- paste(paste('cna.threshold', '<-', yaml$panoply_cmap_analysis$cna_threshold),
                  paste('cna.effects.threshold', '<-', yaml$panoply_cmap_analysis$cna_effects_threshold), 
                  paste('min.sigevents', '<-', yaml$panoply_cmap_analysis$min_sigevents),
                  paste('max.sigevents', '<-', yaml$panoply_cmap_analysis$max_sigevents),
                  paste('top.N', '<-', yaml$panoply_cmap_analysis$top_N),
                  paste('fdr.pvalue', '<-', yaml$panoply_cmap_analysis$fdr_pvalue),
                  paste('log.transform', '<-', yaml$panoply_cmap_analysis$log_transform),
                  paste('must.include.genes', '<-', ifelse(is.null(yaml$panoply_cmap_analysis$must_include_genes), 'NULL', yaml$panoply_cmap_analysis$must_include_genes)),
                  paste('cis.fdr', '<-', yaml$panoply_cmap_analysis$cis_fdr),
                  paste('legacy.score', '<-', yaml$panoply_cmap_analysis$legacy_score),
                  paste('rankpt.n', '<-', yaml$panoply_cmap_analysis$rankpt_n),
                  paste('mean.rankpt.threshold', '<-', yaml$panoply_cmap_analysis$mean_rankpt_threshold),
                  paste('cmap.fdr', '<-', yaml$panoply_cmap_analysis$cmap_fdr),
                  paste('alpha', '<-', yaml$panoply_cmap_analysis$alpha),
                  sep = "\n")
  write(output, 'cmap-config-custom.r')
}

########################################################################################

# Main function and logic
# How to handle the commandline variables based on current module?
parse_command_line_parameters <- function(opt){
  #If there are specific command line values given:
  #Get the module, read the yaml, and parse/overwrite with command line accordingly
  yaml <- read_yaml_file(opt$master_yaml)
  if (opt$module == 'parse_sm_table' & check_if_any_command_line(opt)){ # TRUE if command line changes are being made
    yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed
    yaml <- check_parse_sm_params(opt, yaml) #Returns updated yaml if module params were changed
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
    
  }else if (opt$module == 'normalize_ms_data' & check_if_any_command_line(opt)){
    yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed via command line
    yaml <- check_normalize_sm_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
    
  }else if (opt$module == 'rna_protein_correlation' & check_if_any_command_line(opt)){
    #yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed via command line
    yaml <- check_rna_protein_correlation_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
    
  }else if (opt$module == 'harmonize' & check_if_any_command_line(opt)){
    yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed via command line
    yaml <- check_harmonize_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
    
  }else if (opt$module == 'sample_qc' & check_if_any_command_line(opt)){
    #yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed via command line
    yaml <- check_sample_qc_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
    
  }else if (opt$module == 'cna_analysis' & check_if_any_command_line(opt)){
    #yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed via command line
    yaml <- check_cna_analysis_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
    
  }else if (opt$module == 'cna_corr_report' & check_if_any_command_line(opt)){
    #yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed via command line
    yaml <- check_cna_corr_report_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
   
    
  }else if (opt$module == 'rna_corr_report' & check_if_any_command_line(opt)){
    #yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed via command line
    yaml <- check_rna_corr_report_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
   
    
  }else if (opt$module == 'cons_clust' & check_if_any_command_line(opt)){
    #yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed via command line
    yaml <- check_cons_clust_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
    
  }else if (opt$module == 'association' & check_if_any_command_line(opt)){
    yaml <- check_global_params(opt, yaml) #Returns updated yaml if globals were changed via command line
    yaml <- check_association_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
  }else if (opt$module == 'association_report' & check_if_any_command_line(opt)){
    yaml <- check_association_report_params(opt,yaml)
    write_custom_config(yaml)
    
  }else if (opt$module == 'cmap_analysis'){
    yaml <- check_cmap_analysis_params(opt, yaml) #Returns updated yaml if module params were changed via command line
    write_custom_cmap_config(yaml) #Write params to custom-config.r (CMAP version)
   
    
  }else if (opt$module == 'immune_analysis' & check_if_any_command_line(opt)){
    yaml <- check_immune_analysis_params(opt,yaml)
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
  
  }else if (opt$module == 'blacksheep' & check_if_any_command_line(opt)){
    yaml <- check_blacksheep_params(opt,yaml)
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
  }else if (opt$module == 'mo_nmf' & check_if_any_command_line(opt)){
    yaml <- check_mo_nmf_params(opt,yaml)
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
  }else if (opt$module == 'ssgsea_projection' & check_if_any_command_line(opt)){
    yaml <- check_ssgsea_projection_params(opt,yaml)
    write_custom_config(yaml) #Write params to custom-config.r (GENERIC)
    
  }else if (opt$module == 'pipeline'){
    yaml <- check_pipeline_params(opt, yaml)
    write_custom_config(yaml)
    write_custom_cmap_config(yaml)
    #write out final params for user
  }else{
    #There were no command line values given for any of the given modules
    write_custom_config(yaml) #Write params to custom-config.r  (GENERIC)
    
  }
  write_yaml_file(yaml) #write out final params for user (GENERIC) can be made more specific if needed
}


# MAIN:
parse_command_line_parameters(opt)


