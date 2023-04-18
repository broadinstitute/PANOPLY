#' @title Extract mutations and sample IDs into a long format using transcript IDs
#'
#' @param maf_file Path to .maf file with mutation calls
#' @param transcript_col Column name for chosen mapped transcript
#' @param other_transcript_col Column name for other mapped transcript
#'
#' @return Dataframe  | gene (transcript) | sample_id | position | wt_residue | mut_residue |
get_mutations_all <- function(
    maf_file,
    sample_col = "Sample.ID",
    transcript_col = "Annotation_Transcript",
    other_transcript_col = "Other_Transcripts",
    variant_types = c("SNP"),
    variant_classes = c("Missense_Mutation")
) {
  wxs_maf <- read.maf(maf = maf_file)
  
  cols <- c(sample_col, transcript_col, other_transcript_col, "Variant_Type", "Variant_Classification", "Protein_Change")
  all_mutations <- wxs_maf@data[, mget(cols)]
  filtered_muts <- all_mutations %>% filter(Variant_Type %in% variant_types & Variant_Classification %in% variant_classes)
  
  # TODO: write out examples for mapping
  # grab all other transcripts and corresponding protein changes
  all_other_transcript_mutations <- str_extract_all(filtered_muts$Other_Transcripts, "ENST\\d+.\\d+_\\w*_p.[A-Z]\\d+[A-Z]")
  all_other_transcripts <- str_extract_all(all_other_transcript_mutations, "ENST\\d+.\\d+")
  all_other_protein_change <- str_extract_all(all_other_transcript_mutations, "[A-Z]\\d+[A-Z]")
  protein_change <- str_extract_all(filtered_muts$Protein_Change, "[A-Z]\\d+[A-Z]")
  all_transcripts <- list()
  all_protein_change <- list()
  
  for (i in 1:length(all_other_transcripts)) {
    all_transcripts[i] <- paste(c(all_other_transcripts[[i]], filtered_muts$Annotation_Transcript[i]), collapse = "|")
    all_protein_change[i] <- paste(c(all_other_protein_change[[i]], protein_change[[i]]), collapse = "|")
  }
  filtered_muts$all_transcripts <- unlist(all_transcripts)
  filtered_muts$all_protein_change <- unlist(all_protein_change)
  
  all_mutations_exploded <- filtered_muts %>% separate_rows(all_transcripts, all_protein_change)
  
  mutations <- all_mutations_exploded[, "all_transcripts"]
  colnames(mutations) <- c("gene")
  mutations$sample_id <- all_mutations_exploded[[sample_col]]
  
  mutations$position <- as.numeric(str_extract(all_mutations_exploded$all_protein_change, "(\\d+)"))
  mutations$wt_residue <- str_extract(all_mutations_exploded$all_protein_change, "[A-Z]")
  mutations$mut_residue <- str_extract_all(all_mutations_exploded$all_protein_change, "[A-Z]", simplify = TRUE)[, 2]
  
  return(as.data.frame(mutations))
}


#' @title Obtain all PTM sites from GCT file
#'
#' @param ptm_gct Path to .gct file with PTM
#' @param level Level at which to annotate site locations. Options: "protein" (default), "gene", "transcript"
#'
#' @return Dataframe  | gene (protein) | position | residue | kinase |
get_ptm_sites <- function(ptm_gct, id_column = "ProteinID", residues = "STYK") {
  ptm <- parse_gctx(ptm_gct)
  rdesc <- ptm@rdesc
  
  mod_res_pattern <- paste0("[", residues, "]\\d+")  # regex pattern "[STYK]\\d+"
  variable_sites <- str_extract_all(rownames(rdesc), mod_res_pattern)
  names(variable_sites) <- rdesc[, id_column]
  id_variable_sites <- unlist(variable_sites)
  
  ptm_sites <- as.data.frame(names(id_variable_sites))
  colnames(ptm_sites) <- c("gene")
  
  ptm_sites$position <- as.numeric(str_extract(id_variable_sites, "(\\d+)"))
  ptm_sites$residue <- str_extract(id_variable_sites, "[A-Z]")
  ptm_sites$kinase <- ""  # placeholder to run ActiveDriver
  
  return(ptm_sites)
}


#' @title Format FASTA protein sequences into ActiveDriver format
#'
#' @param proteome_fasta Path to .gct file with PTM
#' @param level Level at which to annotate sequences. Options: "protein" (default), "gene", "transcript"
#'
#' @return ActiveDriver formatted list with level as names and AA sequences as value
get_sequences <- function(proteome_fasta, level = "protein") {
  sequences <- proteome_fasta
  names(sequences) <- extract_names(proteome_fasta, level)
  sequences <- unlist(sequences)
  
  return(sequences)
}

#' @title Format FASTA protein sequences into ActiveDriver format
#'
#' @param sequences ActiveDriver formatted list of sequences
#'
#' @return ActiveDriver formatted list with level as names and zeroes as value
get_disordered_regions_zero <- function(sequences) {
  disordered_sequences <- list()
  
  for (gene in names(sequences)) {
    template <- strrep("0", nchar(sequences[gene]))
    disordered_sequences[gene] <- template
  }
  disordered_sequences <- unlist(disordered_sequences)

  return(disordered_sequences)
}


#' @title mapping from transcript to protein and gene IDs
#'
#' @param proteome_fasta path to FASTA protein sequences
#'
#' @return transcript | protein | gene_id |
get_transcript_map <- function(proteome_fasta) {
  protein <- extract_names(proteome_fasta, "protein")
  transcript <- extract_names(proteome_fasta, "transcript")
  gene <- extract_names(proteome_fasta, "gene")
  geneSymbol <- extract_names(proteome_fasta, "geneSymbol")
  
  transcript_map <- data.frame(protein, gene, geneSymbol)
  rownames(transcript_map) <- transcript
  
  return(transcript_map)
}


#' @title mapping from protein to transcript and gene IDs
#'
#' @param transcript_map
#'
#' @return protein | transcript | gene_id |
get_protein_map <- function(transcript_map) {
  protein_map <- transcript_map
  rownames(protein_map) <- protein_map$protein
  protein_map$transcript <- rownames(transcript_map)
  
  return(protein_map)
}

extract_names <- function(fasta, level = "protein") {
  if (level == "gene") {
    ids_regex <- "ENSG\\d+.\\d+"
  } else if (level == "protein") {
    ids_regex <- "ENSP\\d+.\\d+"
  } else if (level == "transcript") {
    ids_regex <- "ENST\\d+.\\d+"
  } else if (level == "geneSymbol") {
    ids_regex <- "GN=\\w+"
  } else {
    stop('Supplied naming level does not exist. Options: "gene", "protein", "transcript"')
  } 
  seq_ids <- str_extract(names(fasta), ids_regex)
  
  return(seq_ids)
}

AD_extract_top_genes <- function(results, protein_map = NULL, p_thres = 0.05, sort = "fdr") {
  gene_pval <- results$all_gene_based_fdr
  if (!is.null(protein_map)) {
    gene_pval$geneSymbol <- protein_map[gene_pval$gene, "geneSymbol"]
  }
  # gene_pval_filt <- gene_pval %>% filter(p <= p_thres)
  gene_ranked <- gene_pval[order(unlist(gene_pval[, sort])), ]
  
  return(gene_ranked)
}



##### OLD FUNCTIONS #####
#' @title Extract mutations and sample IDs into a long format
#'
#' @param maf_file Path to .maf file with mutation calls
#' @param gene_col Column name to extract gene/protein/transcript name
#'
#' @return Dataframe | gene | sample_id | position | wt_residue | mut_residue |
get_mutations <- function(
    maf_file,
    gene_col = "Hugo_Symbol",  # Hugo_Symbol, HGNC_UniProt_ID_supplied_by_UniProt_
    variant_types = c("SNP"),
    variant_classes = c("Missense_Mutation")
) {
  
  wxs_maf <- read.maf(maf = maf_file)
  plotmafSummary(maf = wxs_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  
  cols <- c("Sample.ID", gene_col, "Variant_Type", "Variant_Classification", "Protein_Change")
  all_mutations <- wxs_maf@data[, mget(cols)]
  filtered_muts <- all_mutations %>% filter(Variant_Type %in% variant_types & Variant_Classification %in% variant_classes)
  
  mutations <- filtered_muts[, mget(gene_col)]
  colnames(mutations) <- c("gene")
  mutations$sample_id <- filtered_muts$Sample.ID
  
  prot_change <- filtered_muts$Protein_Change
  mutations$position <- as.numeric(str_extract(prot_change, "(\\d+)"))
  mutations$wt_residue <- str_extract(prot_change, "[A-Z]")
  mutations$mut_residue <- str_extract_all(prot_change, "[A-Z]", simplify = TRUE)[, 2]
  
  return(mutations)
}


get_ptm_sites_SM_onesite <- function(ptm_gct, level = "protein") {
  ptm <- parse_gctx(ptm_gct)
  if (level == "gene") {
    level_col <- "geneID"
  } else if (level == "protein") {
    level_col <- "id.description"
  } else if (level == "transcript") {
    level_col <- "transcriptID"
  } else {
    stop('Supplied naming level does not exist. Options: "gene", "protein", "transcript"')
  }
  
  rdesc <- ptm@rdesc
  ptm_sites <- as.data.frame(rdesc[, level_col])
  colnames(ptm_sites) <- c("gene")
  
  variable_sites <- rdesc$variableSites
  ptm_sites$position <- as.numeric(str_extract(variable_sites, "(\\d+)"))
  ptm_sites$residue <- str_extract(variable_sites, "[A-Z]")
  ptm_sites$kinase <- ""  # placeholder to run ActiveDriver
  
  return(ptm_sites)
}