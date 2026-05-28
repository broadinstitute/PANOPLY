#!/usr/bin/env Rscript
#
# Copyright (c) 2024 The Broad Institute, Inc. All rights reserved.
#
rm(list=ls())
options( warn = 1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  #### Input Parameters ####
  make_option( c("-m", "--metabolome_gct"), action='store', type='character',  dest='metabolome_gct', help='GCT file containing expression data for the metabolome.'),
  make_option( c("-n", "--meta_id_type"), action='store', type='character', dest='meta_id_type', help='Type of ID type used for metabolites.'),#, default='hmdb_id'),
  make_option( c("-i", "--meta_id_col"), action='store', type='character', dest='meta_id_col', help='Rdesc column in metabolome GCT containing metabolite IDs.'),#, default='hmdb_id'),
  make_option( c("-o", "--ome_gct"), action='store', type='character',  dest='ome_gct', help='GCT file containing expression data for an additional -ome (e.g. proteome, transcriptome, etc.)'),
  make_option( c("-t", "--ome_type"), action='store', type='character',  dest='ome_type', help='Label for the additional ome-type.'),
  make_option( c("-c", "--gene_column"), action='store', type='character', dest='gene_col', help='Column name in rdesc in the GCT that contains gene names.', default='geneSymbol'),
  make_option( c("-d", "--gene_id_type"), action='store', type='character', dest='gene_id_type', help='Type of ID contianed in gene_col', default='SYMBOL'),
  make_option( c("-g", "--groups_file"), action='store', type='character',  dest='groups_file', help='Groups-file, i.e. an annotations file subsetted to annotations of interest.'),
  #### Analysis ####
  make_option( c("-w", "--pathway_db"), action='store', type='character',  dest='pthw_db', help='Pathway database to use. KEGG pathways will be used by default, if missing.'),
  make_option( c("-l", "--max_annot_levels"), action='store', type='numeric', dest='max_annot_levels', help='Maximum number of levels an annotation can have and be considered discrete.'), # default='10'),
  make_option( c("-a", "--anal_type"), action='store', type='character', dest='anal_type', help='Analysis method to use ("ORA" for Overrepresentation Analysis or "QEA" for Quantitative Enrichment Analysis).'), 
  make_option( c("-b", "--pval_comb"), action='store', type='character', dest='pval_comb', help='Method for combining p-values in multiomic enrichment analysis. Options include "query" (combine queries), "pvalu" (unweighted), "pvalo" (overall), or "pvalp" (pathway-level).'), 
  make_option( c("-p", "--pval_signif"), action='store', type='numeric', dest='pval_signif', help='P-value threshold for significant enrichement.'), 
  make_option( c("-k", "--top_n_networks"), action='store', type='numeric', dest='top_n_networks', help='Top N networks to plot per annot subvalue.'), 
  make_option( c("-r", "--impact_metric"), action='store', type='numeric', dest='impact_metric', help="Topological impact metric to be used in plotting ('Impact.BC' for betweenness centrality, 'Impact.CC' for closeness centrality, or 'Impact.DC' for degree centrality)."), 
  make_option( c("--min_overlap"), action='store', type='numeric', dest='min_overlap', help='Minimum number of overlapping features required for pathway enrichment analysis.'), 
  make_option( c("--background_filter"), action='store', type='logical', dest='background_filter', help='Whether to filter pathway entries to only those present in the dataset (TRUE) or use all pathway entries (FALSE). Default is TRUE.'), 
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.'),  # default = 2),
  make_option( c("-f", "--output_directory"), action='store', type='character',  dest='output_dir', help='Directory to output files to.',  default = 'results/'),
  make_option( c("-y", "--yaml"), action="store", dest='yaml_file', type="character", help="Path to .yaml file with parameters."),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)

#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                      #  # for testing arguments
                      #  args = c(
                      #    # '--metabolome_gct',"/opt/input/HMDB_ID_GCTs/ODG-v2_2-metabolomics_log_norm-HMDB_UNIQUE.gct",
                      #    '--metabolome_gct',"/opt/input/ODG-v4-metabolome-all-log2-median-norm-QCfilter-NMFk3core_n56x647.gct",
                      #    # '-n',"hmdb_id",
                      #    '-i',"HMDB.ID",
                      #    '-w',"smpdb_pathway",
                      #    # '-n',"kegg_id",
                      #    # '-i',"KEGG.ID",
                      #   #  '--ome_gct',"opt/input/ODG-v4-proteome-SpectrumMill-ratio-QCfilter-NArm-NMFk3core_n80x12650.gct",
                      #   #  '-t',"proteome",
                      #    #   # '--ome_gct',"opt/input/ODG-v2_2-rnaseq-expression-TPM-protein-coding-log2-median-norm-NArm-with-NMF.gct",
                      #    #   # '-t',"RNA",
                      #    #   # '-g',"opt/input/sample-info.csv",
                      #    '-g',"opt/input/groups-subset-subsetted.csv",
                      #    #   '-l',"15",
                      #      '-a',"QEA",
                      #    # '-b',"pvalo",
                      #    #   '-p',"0.01",
                      #    #   '-k',"15",
                      #    # '-r',"Impact.CC",
                      #    # '--background_filter', 'FALSE',
                      #    # '--min_overlap', '0',
                      #    '-y',"opt/input/master-parameters.yaml",
                      #    #   # '-f',"/opt/input/prelim_results",
                      #    '-x',"ODG_v3")
)

#### Parse YAML Arguments ####
opt = opt_cmd # initialize options with command line options
if ( !is.null(opt$yaml_file) ) {
  #### read in yaml ####
  library(yaml)
  yaml_out <- read_yaml(opt$yaml_file)
  #### locate the NMF-postprocessing parameters section ####
  if ( !is.null(yaml_out$panoply_metaboanalyst) ) { # if we have an NMF parameters section
    yaml_metaboanalyst =  yaml_out$panoply_metaboanalyst # read in those parameters
  } else { # otherwise, stop
    stop(glue("The parameter file '{opt$yaml_file}' does not contain a MetaboAnalyst parameters section. Please check that the yaml file contains the appropraite 'panoply_metaboanalyst' section."))
  }
  #### overwrite the command-line parameters ####
  # global parameters
  if (is.null(opt$gene_col)) opt$gene_col = yaml_out$global_parameters$gene_mapping$gene_id_col
  # postprocessing parameters
  if (is.null(opt$meta_id_type)) opt$meta_id_type = yaml_metaboanalyst$meta_id_type
  if (is.null(opt$meta_id_col)) opt$meta_id_col = yaml_metaboanalyst$meta_id_col
  if (is.null(opt$gene_id_type)) opt$gene_id_type = yaml_metaboanalyst$gene_id_type
  if (is.null(opt$max_annot_levels)) opt$max_annot_levels = yaml_metaboanalyst$max_annot_levels
  if (is.null(opt$anal_type)) opt$anal_type = yaml_metaboanalyst$anal_type
  if (is.null(opt$pval_comb)) opt$pval_comb = yaml_metaboanalyst$pval_comb
  if (is.null(opt$pval_signif)) opt$pval_signif = yaml_metaboanalyst$pval_signif
  if (is.null(opt$top_n_networks)) opt$top_n_networks = yaml_metaboanalyst$top_n_networks
  if (is.null(opt$impact_metric)) opt$impact_metric = yaml_metaboanalyst$impact_metric
  if (is.null(opt$pthw_db)) opt$pthw_db = yaml_metaboanalyst$pthw_db
  if (is.null(opt$min_overlap)) opt$min_overlap = yaml_metaboanalyst$min_overlap
  if (is.null(opt$background_filter)) opt$background_filter = yaml_metaboanalyst$background_filter
} else { # if no YAML was provieded
  # check if any necessary parameters are missing
  if( any(sapply(list(opt$meta_id_type, 
                      # opt$gene_col, opt$gene_id_type, # allow gene_col and gene_id_type to be missing, since omic GCT is optional
                      opt$max_annot_levels, opt$anal_type, opt$pval_comb, opt$pval_signif, opt$top_n_networks), 
                 is.null)) ) { # if we have at least one missing parameter
    stop("Master Parameter yaml-file is missing. Please either provide a master-parameters file, or manually provide all NMF parameters.") # error and stop
  }
}
# ensure output directory exists
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir) # create output directory if it does not exist
# save opt to .Rdata object
save(file = file.path(opt$output_dir, "metaboanalyst_opt.Rdata"), opt)


library(glue)
library(cmapR)
library(tidyverse)
library(globaltest)
library(ActivePathways)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(igraph)
library(d3r)
library(readxl)
library(ComplexHeatmap)

library(rjson)
library(jsonlite)

library(SimDesign) # for quiet() function

overwrite_rid = function(gct, new_rid, allow_dups=F) {
  if (length(new_rid) != length(gct@rid)) stop("RID do not match in length")
  if (!allow_dups && length(new_rid) != length(unique(new_rid))) stop("RID are not unique")
  gct@rid = new_rid
  gct@rdesc$id = new_rid
  rownames(gct@mat) = new_rid
  return(gct)
}


max_hm_pathways = 50
print_internal_placemarks = FALSE # toggle to print placemarkers for Enrichment Analysis on each Annotation Subvalue (will not suppress warnings)

################################
####      Read In Data      ####
################################
cat("\n\n####################\nReading In Data\n\n")
#### Read in Omics Data ####
gct_meta = parse_gctx(opt$metabolome_gct)
# optionally read in additional ome
if (!is.null(opt$ome_gct)) {
  multiomic = TRUE
  gct_ome = parse_gctx(opt$ome_gct)
} else {
  multiomic = FALSE
  gct_ome = NULL
}


#### Read in Annotations & ID Mapping ####
if (!is.null(opt$groups_file)) { annots = read.csv(opt$groups_file) } else { # read in groups file or
  annots = gct_meta@cdesc; cat("\nWARNING: No groups file provided; module will attempt to use metabolome gct@cdesc.\n") } # print warning and use cdescs

compound_map = qs::qread(file.path(opt$lib_dir, "pathway_db/master_compound_db.qs")) # MetaboAnalyst master compound table (ID columns + `name`)
## create stable row IDs
if (!"name" %in% names(compound_map)) stop("master_compound_db is missing required column `name`.")
rownames(compound_map) = make.names(compound_map$name, unique = TRUE) # set rownames
## print warning for duplicate IDs
dup_row_idx = duplicated(compound_map$name) | duplicated(compound_map$name, fromLast = TRUE) # identify duplicated rows
if (any(dup_row_idx)) { # print warning about duplicated metabolite names
  example = paste(head(rownames(compound_map)[dup_row_idx]), collapse = ";\n")
  warning(glue("master_compound_db: {sum(dup_row_idx)} row(s) share a duplicated `name`; `master_id` from make.names(..., unique=TRUE) disambiguates them (e.g. suffix .1, .2).\nExample name(s):\n{example}"))
}

if (is.null(opt$pthw_db)) {
  warning("No pathway database selected. Defaulting to KEGG database.")
  opt$pthw_db = "kegg_pathway"
}
kegg_dbs = c('kegg_pathway'='hsa', 'kegg_pathway_2023'='hsa_2023') # kegg db files
valid_qs_files = c(names(kegg_dbs), # manually add Kegg DBs pathway names; these are handled manually since they require different wrangling
                   setdiff(gsub('.qs$','',list.files( file.path(opt$lib_dir, "pathway_db"), pattern = '.qs')),
                           c(kegg_dbs, names(kegg_dbs), # exclude the multiomic KEGG DB *files*; these are handled manually since they require different wrangling
                             c('compound_db', 'lipid_compound_db', 'master_compound_db')))) # exclude compound database files; these aren't pathways
if (!opt$pthw_db %in% valid_qs_files) {
  stop(glue("Invalid pathway database '{opt$pthw_db}'. Please select one of the following:\n{paste0(valid_qs_files, collapse='\n' )}"))
}
# do not run multiomic analysis with the SMP Database
if (multiomic && !( opt$pthw_db %in% kegg_dbs ) ) {
  stop(glue("This pathway ('{opt$pthw_db}') contains only metabolites, and does not support multiomic analysis. Please use a KEGG database ('{paste(names(kegg_dbs), collapse='\\' or \\'')}') instead."))
}

#### Read in / Format Pathways ####
if (opt$pthw_db %in% names(kegg_dbs)) { # if we're looking at the multiomic KEGG databases
  # format QS object from MetaboAnalyst Site into GMT-adjacent
  pathways_qs = qs::qread(file.path(opt$lib_dir,"pathway_db", paste0(kegg_dbs[opt$pthw_db], ".qs")))
  pathways = lapply(unname(pathways_qs$path.ids), function(pathway) {
    pathway_list = list(ID = pathway,
                        name = names(pathways_qs$path.ids[pathways_qs$path.ids==pathway]),
                        entries = unlist(sapply(pathways_qs$mset.list[[pathway]], USE.NAMES=F, # unlist pathway elements
                                                function(x) {strsplit(x, " ", fixed = TRUE)} )))#, # and separate strings into single entries
                        # cmpd.counts = sum(grepl(glue("^(cpd)|(gl):"), pathways_qs$mset.list[[pathway]])), # count compounds (before unlisting)
                        # gene.counts = sum(grepl(glue("^hsa:"), pathways_qs$mset.list[[pathway]]))) # count genes (before unlisting)
    return(pathway_list)
  })
  names(pathways) = pathways_qs$path.ids
  # set pathway IDs
  pathway_id_type = list(meta="kegg_id",
                         gene="ENTREZID")
  
  # # get unique compound / gene counts
  # unique.cmpd = pathways_qs$uniq.cmpd.count
  # unique.gene = pathways_qs$uniq.gene.count
  
  # binary toggle for topological measures
  has_topology = TRUE
  pathway_topology_scores = lapply(names(pathways), function(pathway) {
    pathway_list = list(ID = pathway,
                        entries.grouped = pathways_qs$mset.list[[pathway]],
                        bc.list = pathways_qs$bc.list[[pathway]],
                        dc.list = pathways_qs$dc.list[[pathway]],
                        cc.list = pathways_qs$cc.list[[pathway]])
  })
  names(pathway_topology_scores) = names(pathways)
  
  # graph list
  graph_list = pathways_qs$graph.list
} else {
  # # for loop is JUST for testing databases; should be commented out for general usage 
  # for (pthw_db in setdiff(valid_qs_files, c(names(kegg_dbs)))) { # check that all QS files work with this format
  #   opt$pthw_db = pthw_db
  #   cat("\n\n###################\n\n")
  #   cat(pthw_db); cat('\n\n')
    
    #   wrangle .qs dataframes into a list structure
    pathways_qs = qs::qread(file.path(opt$lib_dir,"pathway_db", paste0(opt$pthw_db, ".qs")))
    pathways = apply(pathways_qs, 1,
      function(p) {
        list(ID = p['id'][[1]], #p['image'][[1]], # image is only valid for a few DBs
            name = p['name'][[1]],
            reference = p['reference'][[1]],
            entries = str_split(p['member'], "; ")[[1]])
      })
    names(pathways) = pathways_qs$id #pathways_qs$image
    # drop malformed pathways
    pthw_keep = names(which(!unlist(lapply(pathways, function(p) {any(is.na(p$entries))})))) # prune pathways with NA features
    if ( length(pthw_keep) != length(pathways) ) {
      pthw_drop = setdiff(names(pathways), pthw_keep)
      warning(glue("Dropping pathway(s) from '{opt$pthw_db}' database with NA entries:\n{paste(pthw_drop, paste0(tidyr::replace_na(pathways[[pthw_drop]]$entries, 'NA'), collapse=', '), sep='; Entries: ', collapse='\n')}"))
    }
    pathways = pathways[pthw_keep]

    # set pathway IDs
    pathway_id_type = list(meta="name")
    pthw_entries = unique(unlist(lapply(pathways, `[[`, "entries")))
    if (any(!pthw_entries %in% compound_map[[pathway_id_type$meta]])) {
      warning(glue("Pathway entries missing from master_compound_db `{pathway_id_type$meta}`:\n{paste0(pthw_entries[which(!pthw_entries %in% compound_map[[pathway_id_type$meta]])], collapse=', ')}"))
    }
  #   cat("\n\n###################\n\n")
  #   # print(head(pathways_qs))
  #   # print(head(pathways_qs$id))
  #   # print(head(pathways_qs$image))
  #   print(pathways[[1]])
  # }

  # pathways_gmt = read.GMT(file.path(opt$lib_dir, "pathway_db/smpdb_pathway.gmt"))
  # pathways = lapply(pathways_gmt, function(p) {
  #   pathway_list = list(ID = p$id,
  #                       name = p$name,
  #                       entries = p$genes)#,
  #                       # cmpd.counts = length(unique(p$genes)) # only take unique values
  #                       #               + max(sum(p$genes=='NA')-1, 0) ) # but also count every NA as a unique compound
  #   return(pathway_list)
  # })
  
  # # get unique compound counts
  # unique.cmpd = length(unique(unlist(lapply(pathways_gmt, function(p) {p$genes})))) # get unique genes
  
  # binary toggle for topological measures
  has_topology = FALSE
  # graph list
  graph_list = NULL
}

#### Initialize Network Annotation list ####
network_annots = list()


################################
#### Map ID / Pathway Match ####
################################


cat("\n\n####################\nMetabolite ID Validation\n\n")
if (! opt$meta_id_type %in% names(compound_map)) stop(glue("Provided metabolic ID type '{opt$meta_id_type}' is not one of the accepted ID types ({paste(names(compound_map), collapse=', ')})"))
if (! pathway_id_type$meta %in% names(compound_map)) stop(glue("Mapped metabolic ID type '{pathway_id_type$meta}' is not one of the accepted ID types ({paste(names(compound_map), collapse=', ')})"))

# get relevant IDs from rid or rdesc
if (is.null(opt$meta_id_col)) { cpd_vec = gct_meta@rid } else { cpd_vec = gct_meta@rdesc[[opt$meta_id_col]] }

valid_cpd_rid = gct_meta@rid[which( !is.na(cpd_vec) & (cpd_vec %in% compound_map[[opt$meta_id_type]]) )]
if (length(valid_cpd_rid)==0) stop(glue("No IDs in the GCT mapped to valid compounds. Please check that your data uses {opt$meta_id_type} IDs, or select a different ID type."))
# subset to valid compound IDs
# toDo: add lipid ID mapping
cat(glue("\nOut of {length(cpd_vec)} features, {length(valid_cpd_rid)} appear in our DB as valid {opt$meta_id_type}.\n\n"))
gct_meta_filt = subset_gct(gct_meta, rid=valid_cpd_rid)
# map to pathway ID type (resolve duplicate lookup rows: first non-NA target, else drop feature)
if (opt$meta_id_type != pathway_id_type$meta || # if we need to change ID type
    !is.null(opt$meta_id_col)) { # OR the id column wasn't the rid
  
  # helper function-- if compound_map has multiple rows per input ID, use the first non-NA pathway target ID.
  resolve_pathway_meta_ids = function(input_ids, compound_map, from_col, to_col) {
    lookup = compound_map[[from_col]]
    target = compound_map[[to_col]]
    out = vapply(input_ids, function(id) {
      if (is.na(id) || !nzchar(trimws(as.character(id)))) return(NA_character_)  # if our ID is NA, return NA
      hits = which(!is.na(lookup) & lookup == id) # locate all plausible ID matches
      if (length(hits) == 0L) return(NA_character_) # if we have no hits, return NA
      vals = target[hits] # get hits
      ok = !is.na(vals) & nzchar(trimws(as.character(vals))) # check for a valid hit
      if (!any(ok)) return(NA_character_) # if we have no valid hits, return NA
      as.character(vals[which(ok)[1L]]) # get first valid hit
    }, character(1), USE.NAMES = FALSE)
    u_ids = unique(input_ids[!is.na(input_ids) & nzchar(trimws(as.character(input_ids)))])
    n_multi = sum(vapply(u_ids, function(id) sum(!is.na(lookup) & lookup == id, na.rm = TRUE) > 1L, logical(1)))
    if (n_multi > 0L) {
      warning(glue("{n_multi} input `{from_col}` ID(s) has multiple entries in the metabolite mapping database; using the first non-NA `{to_col}` per ID."))
    }
    out
  }

  # get relevant IDs from rid or rdesc
  if (is.null(opt$meta_id_col)) { cpd_vec_filt = gct_meta_filt@rid } else { cpd_vec_filt = gct_meta_filt@rdesc[[opt$meta_id_col]] }
  new_rid = resolve_pathway_meta_ids(cpd_vec_filt, compound_map, opt$meta_id_type, pathway_id_type$meta)
  has_pathway_id = !is.na(new_rid) & nzchar(new_rid)
  if (any(!has_pathway_id)) {
    warning(glue("Excluding {sum(!has_pathway_id)} feature(s): matched on `{opt$meta_id_type}` but no non-NA `{pathway_id_type$meta}` ID."))
    # gct_meta_filt = subset_gct(gct_meta_filt, rid = gct_meta_filt@rid[which(has_pathway_id)])
    # new_rid = new_rid[has_pathway_id]
    # cpd_vec_fin = cpd_vec_filt[has_pathway_id]
  }
  rid_dup = !(new_rid %in% unique(new_rid[duplicated(new_rid)])) # identify duplicated RIDs (covers NA values)
  # overwrite RID and drop duplicated
  tmp = overwrite_rid(gct_meta_filt, new_rid, allow_dups = T)
  meta_val = subset_gct(tmp, rid_dup)
  cat(glue("\nOut of {length(cpd_vec_filt)} features, {length(meta_val@rid)} had valid and unique {pathway_id_type$meta} IDs.\n\n"))
} else { meta_val = gct_meta_filt }
# initialize feature_map map between pathway IDs and human-readable names (used throughout)
feature_map = compound_map$name # get compound names
names(feature_map) = compound_map[, c(pathway_id_type$meta)] # name with pathway IDs
# optionally add KEGG formatting (cpd: prefix) to metabolite IDs
if (opt$pthw_db %in% names(kegg_dbs) || multiomic) { # if KEGG DB is used, or if multiomic analysis is performed
  meta_val = overwrite_rid(meta_val, paste0("cpd:", meta_val@rid))
  names(feature_map) = paste0("cpd:", names(feature_map)) # add cpd prefix to feature_map
}

if (multiomic) {
  cat("\n\n####################\nGene ID Validation\n\n")
  keytypes = AnnotationDbi::keytypes(org.Hs.eg.db)
  if (! opt$gene_id_type %in% keytypes) stop(glue("Provided ID type '{opt$gene_id_type}' is not one of the accepted ID types ({paste(keytypes, collapse=', ')})"))
  if (! pathway_id_type$gene %in% keytypes) stop(glue("Mapped ID type '{pathway_id_type$gene}' is not one of the accepted ID types ({paste(keytypes, collapse=', ')})"))
  gene_map = AnnotationDbi::select(org.Hs.eg.db, gct_ome@rdesc[[opt$gene_col]], pathway_id_type$gene, opt$gene_id_type)
  if (dim(gene_map)[1]==0) stop(glue("No {opt$ome_type} features had valid {opt$gene_id_type} IDs in the rdesc '{opt$gene_col}' column."))
  # map to new ID type if necessary
  if (opt$gene_id_type != pathway_id_type$gene) {
    new_rid = gene_map[match(gct_ome@rdesc[[opt$gene_col]], gene_map[[opt$gene_id_type]]), # match GCT rid to compound map
                       pathway_id_type$gene] # and select new ID type
    rid_dup = !(new_rid %in% unique(new_rid[duplicated(new_rid)])) # identify duplicated RIDs (covers NA values)
    # overwrite RID and drop duplicated
    tmp = overwrite_rid(gct_ome, new_rid, allow_dups = T) 
    ome_val = subset_gct(tmp, rid_dup)
    cat(glue("\nOut of {length(gct_ome@rid)} features, {length(ome_val@rid)} had valid and unique {pathway_id_type$gene} IDs.\n\n"))
  } else { ome_val = gct_ome }
  
  cat("\n\n####################\nAdding Multi-Omic Formatting to IDs\n\n")
  meta_mo = meta_val # IDs have already been formatted, if KEGG DB is used
  ome_mo = overwrite_rid(ome_val, paste0("hsa:", ome_val@rid))
  
  # map between IDs and human-readable names
  gene_sym_map_df = AnnotationDbi::select(org.Hs.eg.db, ome_val@rid, 'SYMBOL', pathway_id_type$gene) # get map between pathway IDs and geneSymbols
  gene_sym_map = gene_sym_map_df$SYMBOL # get geneSymbols
  names(gene_sym_map) = paste0("hsa:", gene_sym_map_df[,pathway_id_type$gene]) # match with pathway IDs (with hsa: prefix)
  # append gene IDs to feature_map
  feature_map = c(feature_map, gene_sym_map)

}
# append all pathway gene IDs to feature_map, for use in network diagrams
if (opt$pthw_db %in% names(kegg_dbs)) {
  pathway_gene_ids = unique(unlist(lapply(pathways, function(p) {
    p$entries[grepl("^hsa:", p$entries)]
  }))) %>% gsub("^hsa:", "", .)
  pathway_gene_sym_map_df = AnnotationDbi::select(org.Hs.eg.db, pathway_gene_ids, 'SYMBOL', pathway_id_type$gene)
  pathway_feature_map = pathway_gene_sym_map_df$SYMBOL
  names(pathway_feature_map) = paste0("hsa:", pathway_gene_sym_map_df[, pathway_id_type$gene])
  feature_map = c(feature_map, pathway_feature_map)
  feature_map = feature_map[!duplicated(names(feature_map))] # drop duplicated features
}

#### Write Networks to JSON ####
if (!is.null(graph_list)) {
  pthw_json = list()
  for (pthw in names(graph_list)) {
    pthw_json[[pthw]] = d3r::d3_igraph(igrf = graph_list[[pthw]], json = FALSE)

    # Convert node names to readable labels for network diagram
    node_entries = sapply(pthw_json[[pthw]]$nodes$names, function(node_name) {
      strsplit(node_name, " ", fixed = TRUE) %>% unlist()
    }) # get raw IDs from node names
    pthw_json[[pthw]]$nodes$graphics_name = sapply(node_entries, function(entries) {
      entries_hr = unname(feature_map[entries]) %>% as.character()
      missing_idx = (is.na(entries_hr) | entries_hr == "" ) # identify indices with missing IDs
      entries_hr[missing_idx] = entries[missing_idx] # use original ID if no human-readable name is available
      paste(unique(entries_hr), collapse = ", ") # collapse unique IDs into a single string
    }) # build readable node labels from all compounds/genes in each node
  }
  
  exportJSON <- toJSON(pthw_json)
  write_json(exportJSON, path="networks.JSON")
  file.copy('networks.JSON', opt$output_dir) # copy network JSON to output directory
}



# choose datasets for inputs
if (multiomic) {
  meta_input = meta_mo
  ome_input = ome_mo
} else {
  meta_input = meta_val
}

#### get full universe of IDs ####
all_features = meta_input@rid
if (multiomic) { all_features = c(all_features, ome_input@rid) }



# set up background filtering so that pathways are only tested on features present in the dataset
if ( !opt$background_filter ) {
  # legacy behavior: subset pathways to compound/gene type, but use all IDs from pathway database (not filtered to dataset)
  
  ## get feature counts for all pathways
  if (multiomic) {
    # extract all unique compound and gene IDs from pathway database
    all_pathway_compounds = unique(unlist(lapply(pathways, function(p) {
      p$entries[grepl("^(cpd|gl):", p$entries)]
    })))
    all_pathway_genes = unique(unlist(lapply(pathways, function(p) {
      p$entries[grepl("^hsa:", p$entries)]
    })))
    all_pathway_features = unique(unlist(lapply(pathways, function(p) {p$entries})))
  } else {
    # in single-omic mode, all pathway entries are compounds
    all_pathway_compounds = unique(unlist(lapply(pathways, function(p) {p$entries})))
    all_pathway_genes = NULL
    all_pathway_features = all_pathway_compounds
  }

  # set up background filtering to assume all features in pathway DB are present
  if (multiomic) {
    background_meta = all_pathway_compounds # all compounds in pathway database
    background_ome = all_pathway_genes # all genes in pathway database
    background_mo = all_pathway_features # all features in pathway database
  } else {
    background_meta = all_pathway_compounds # all compounds in pathway database
    background_ome = NULL
    background_mo = NULL
  }
} else {
  # default behavior: filter to dataset features AND subset to compound/gene type
  if (multiomic) {
    # in multiomic mode, separate metabolite and gene features from dataset
    background_meta = all_features[grepl("^(cpd|gl):", all_features)] # metabolite features (cpd: or gl:) in dataset
    background_ome = all_features[grepl("^hsa:", all_features)] # gene features (hsa:) in dataset
    background_mo = all_features # all features in dataset for combined analysis
  } else {
    # in single-omic mode, all_features contains only metabolites
    background_meta = all_features
    background_ome = NULL
    background_mo = NULL
  }
}

################################
####  Perform Enrichement   ####
################################

#### Quantitative Enrichment Analysis Function ####
# takes GCT & annot of interest
# returns transposed matrix & binary numeric cls vector
gct.to.qea.input = function(gct, annot_of_interest, value_of_interest, annots = NULL,
                            rm_na = TRUE,
                            write_to_file = F, prefix="results") {
  #### Data Wrangling ####
  if (is.null(annots)) { annots_sorted = gct@cdesc } else {
    if (sum(gct@cid %in% annots[['Sample.ID']]) == 0 ) stop("No samples had matching IDs in provided annotation table.")  # check that we have overlapping IDs
    annots_sorted = annots[match(gct@cid, annots[['Sample.ID']]),] # reorder to match gct@cid
  }
  
  # drop NA samples, if applicable
  if (rm_na) {
    idx.keep = which(! sapply(annots_sorted[[annot_of_interest]], skip.annot)) # exclude any samples that are of a skippable annotation (e.g. NA, "")
    annots_fin = annots_sorted[idx.keep,]
    gct_fin = subset_gct(gct, cid = idx.keep)
  } else {
    annots_fin = annots_sorted
    gct_fin = gct
  }
  
  # create binary CLS
  cls = annots_fin[[annot_of_interest]] %>%
    { ifelse(!is.na(.) & .==value_of_interest, 1,0) } # NA -> no_<value_of_interest>
  
  # pull out matrix and transpose so columns are features
  mat = as.matrix(t(gct_fin@mat))
  # optionally save input dataset to a file
  if (write_to_file) {
    write.csv(data.frame(cls = cls, mat),
              file=glue("{prefix}_QEA_input_{make.names(annot_of_interest)}_{value_of_interest.name}_vs_{value_of_interest.name}_not.csv"))
  }
  return(list(mat = mat,
              cls = cls))
}
q.ea = function(mat, cls, pathways, background_filter=TRUE, p.val.min=2.3233E-11, min_overlap=NULL) {
  if (dim(mat)[1] != length(cls)) stop(glue("Matrix ({dim(mat)[1]}) and CLS vector ({length(cls)}) have differing dimensions. Please ensure that your matrix has samples for rows and columns for features, and that your CLS has an entry for every sample."))
  if (!is.numeric(cls)) stop(glue("CLS vector is not numeric! Please supply a binary numeric vector to annotate your samples."))
  
  #### Enrichment Analysis ####
  
  # number of hits in pathways
  pathway.hits = sapply(pathways, function(p) {
    x = p$entries
    unique(x[x %in% colnames(mat)])
  })
  if ( max(sapply(pathway.hits, length))==0 ) stop("No hits found in any pathway.")

  # calculate universe size (for use in pval combination)
  if (background_filter) {
    universe.size <- unique(unlist(pathway.hits, function(p) {p$entries})) %>%
      length() # get universe size based on hits present in dataset
  } else {
    universe.size = unique(unlist(lapply(pathways, function(p) {p$entries}))) %>%
      length() # get universe size based on all elements in original pathway database
  }
  
  # calculate enrichment
  gt.obj <- globaltest::gt(cls, mat, subsets = pathway.hits)
  gt.res <- globaltest::result(gt.obj)
  
  # wrangle outputs
  res.df = data.frame(ID = names(pathways),
                      Pathway.Name = unname(sapply(pathways, USE.NAMES=F, function(p) {p$name})), # pull name into DF
                      Observed.Entries = sapply(pathway.hits, length)) %>%
    dplyr::mutate(Observed.Stat = gt.res[, 2],
                  Expected.Stat = gt.res[, 3],
                  Enrichment.Ratio = Observed.Stat/Expected.Stat) %>% # Enrichment.Ratio
    dplyr::mutate(Test.Type = "QEA",
                  Raw.P.Value =  gt.res[, 1] %>% ifelse(.==0, p.val.min, .), # add a minimum p-value of 2.3233E-11, to avoid -log(P.Value)
                  neg.Log.P.Value = -log10(Raw.P.Value)) %>%
    dplyr::filter(Observed.Entries>0) %>% # filter to valid hits
    { if (!is.null(min_overlap)) dplyr::filter(., Observed.Entries >= min_overlap) else . } %>% # filter by minimum overlap if specified
    dplyr::mutate(Holm.P.Value = p.adjust(Raw.P.Value, "holm"),
                  BH.P.Value = p.adjust(Raw.P.Value, "fdr")) %>%
    dplyr::arrange(Raw.P.Value) %>% # sort in order of significance
    # dplyr::mutate_if(is.numeric, signif, 5) %>% # filter to 5 sigfigs # DONT round before combination
    column_to_rownames("ID") # add rownames back
  
  return(list(df = res.df,
              hits = pathway.hits,
              universe.size = universe.size))                  
}


#### Hypergeometric Function ####
# takes GCT & annot of interest
# returns DF with significant features and LogFC values
gct.to.ora.input = function(gct, annot_of_interest, value_of_interest, annots = NULL,
                            significance_cutoff = 0.05, rm_na = TRUE,
                            write_to_file = F, prefix="results") {
  
  #### Data Wrangling ####
  if (is.null(annots)) { annots_sorted = gct@cdesc } else {
    if (sum(gct@cid %in% annots[['Sample.ID']]) == 0 ) stop("No samples had matching IDs in provided annotation table.")  # check that we have overlapping IDs
    annots_sorted = annots[match(gct@cid, annots[['Sample.ID']]),] # reorder to match gct@cid
  }
  
  # drop NA samples, if applicable
  if (rm_na) {
    idx.keep = which(! sapply(annots_sorted[[annot_of_interest]], skip.annot)) # exclude any samples that are of a skippable annotation (e.g. NA, "")
    annots_fin = annots_sorted[idx.keep,]
    gct_fin = subset_gct(gct, cid = idx.keep)
  } else {
    annots_fin = annots_sorted
    gct_fin = gct
  }
  
  # create binary CLS
  cls = annots_fin[[annot_of_interest]] %>%
    { ifelse(!is.na(.) & .==value_of_interest, ., glue("{value_of_interest}_not")) } %>% # NA -> no_<value_of_interest>
    factor(levels = c(glue("{value_of_interest}_not"), value_of_interest)) # factor cls with _not first, to force order of comparison

  #### T-Test ####
  require(limma)
  source('/prot/proteomics/Projects/Protigy/modT.R') # sourcing code from Protigy for modT.test.2class(); file is downloaded in panoply_utils docker
  d = rownames_to_column(as.data.frame(gct_fin@mat), "feature_id")
  out = quiet(modT.test.2class(d, 'tmp', groups=cls, id.col = "feature_id")) # wrapped in quiet() to suppress repeated printouts
  out_df = out$output %>%
    dplyr::arrange(P.Value) %>%
    # { if (!is.null(n_feat)) {top_n(., n=n_feat, wt = -P.Value)} else {.} }
    { if (!is.null(significance_cutoff)) { dplyr::filter(., adj.P.Val < significance_cutoff) } else {.} } # filter to significant logFC
  df = data.frame(id = rownames(out_df['logFC']),
                  logFC = round(out_df[['logFC']], 5),
                  P.Value = out_df[['P.Value']],
                  adj.P.Val = out_df[['adj.P.Val']])
  if (write_to_file) {
    write.table(df, file=glue("{prefix}_logFC_signFeat_{make.names(annot_of_interest)}_{value_of_interest.name}_vs_{value_of_interest.name}_not.tsv"),
                sep = '\t', quote=FALSE, row.names = FALSE)
  }
  return(df)
}
o.ea = function(queries, pathways, p.val.min=2.3233E-11, background=NULL, min_overlap=NULL) {
  # calculate original pathway sizes (before filtering)
  pathway.size.full.vec = sapply(pathways, function(p) {length(unique(p$entries))})
  
  # filter pathway entries to background if provided
  if (!is.null(background)) {
    # filter each pathway's entries to only those in background
    pathways_filtered = lapply(pathways, function(p) {
      p$entries_filtered = p$entries[p$entries %in% background]
      return(p)
    })
    # get all unique entries in pathways that are in background
    current.universe <- unique(unlist(lapply(pathways_filtered, function(p) {p$entries_filtered})))
    # recalculate pathway sizes based on filtered entries
    pathway.size.vec = sapply(pathways_filtered, function(p) {length(unique(p$entries_filtered))})
    # recalculate universe size based on filtered entries
    universe.size = length(current.universe)
  } else {
    # use original pathway entries
    current.universe <- unique(unlist(lapply(pathways, function(p) {p$entries}))) # get all unique entries in pathway
    # if we aren't limiting the feature-space (i.e. only hsa: or only cpd: / gl:), calculate size of featurespace
    pathway.size.vec = sapply(pathways, function(p) {sum(length(unique(p$entries)))})
    universe.size = sum(length(current.universe))
  }
  
  # calculate hits
  queries.subset <- queries[queries %in% current.universe] # subset to queries in current.universe
  query.size = length(queries.subset)
  pathway.hits <- lapply(pathways, function(p) { # for each pathway
    y = queries.subset %in% unlist(p$entries) # identify which queries are hits
    queries.subset[y] # return those hits
  })
  pathway.hits.num = sapply(pathway.hits, length) # count the number of hits per pathway
  
  # perform hypergeometric test (equal end of fisher exact test)
  p.val = phyper(pathway.hits.num - 1, pathway.size.vec, universe.size - pathway.size.vec,
                 query.size, lower.tail = F) %>%
    ifelse(.==0, p.val.min, .) # add a minimum p-value (2.3233E-11), to avoid -log(P.Value)=Inf
  
  res.df = data.frame(ID = sapply(pathways, function(p) {p$ID}),
                      Pathway.Name = sapply(pathways, function(p) {p$name}),
                      Observed.Entries = pathway.size.vec, # features in pathway that appear in dataset (filtered to background)
                      N.Hits = pathway.hits.num) %>%
    dplyr::mutate(Expected.Hits = query.size * (pathway.size.vec/universe.size),
                  Enrichment.Ratio = N.Hits/Expected.Hits) %>% # Enrichment.Ratio
    dplyr::mutate(Test.Type = "ORA",
                  Raw.P.Value = p.val,
                  neg.Log.P.Value = -log10(Raw.P.Value)) %>%
    dplyr::filter(N.Hits>0) %>% # filter to valid hits
    { if (!is.null(min_overlap)) dplyr::filter(., Observed.Entries >= min_overlap) else . } %>% # filter by minimum overlap if specified
    dplyr::mutate(Holm.P.Value = p.adjust(Raw.P.Value, "holm"),
                  BH.P.Value = p.adjust(Raw.P.Value, "fdr")) %>%
    dplyr::arrange(Raw.P.Value) %>% # sort in order of significance
    # dplyr::mutate_if(is.numeric, signif, 5) %>% # filter to 5 sigfigs # DONT round before combination
    column_to_rownames('ID')
  
  return(list(df = res.df,
              hits = pathway.hits,
              query.size = query.size,
              universe.size = universe.size))
}

#### Figure Generation Function ####
plot.network = function(g, hits, logFC_df,
                        title,
                        write_to_file = T, prefix="results", out.dir = '.') {
  g = igraph::upgrade_graph(g) # update graph idk
  
  #### Add logFC color-coding info to Graph ####
  # filter logFC dataframe to relevant hits
  logFC_hits = dplyr::filter(logFC_df, id %in% hits)
  # get vertex info as dataframe
  vertex_df = as_data_frame(g, what = "vertices")
  
  # boolean vector of grouped-entries with hits
  logFC_absmax <- apply(vertex_df, 1, function(v, h, l) {
    v_hits = h[stringr::str_detect(v['names'], h)] # identify any hits for THIS vertex
    
    if(length(v_hits)==0) { # if we have no hits for this vertex
      return( data.frame(hit.absmax = NA, logFC=NA)) # return an empty dataframe
    } else { # if we do have hits for this vertex
      out_df = dplyr::filter(l, id %in% v_hits) %>% # filter the logFC dataframe to those hits
        dplyr::rename('hit.absmax'=id) # rename the ID column
      out_df = out_df[which.max(abs(out_df$logFC)), ] # and filter to the hit with the highest logFC value
      return(out_df)
    }
  }, h = hits, l = logFC_hits) %>% do.call(rbind, .) # return a dataframe
  
  V(g)$hit.absmax = logFC_absmax[V(g),'hit.absmax']
  V(g)$logFC = logFC_absmax[V(g),'logFC']
  
  # hits_df = vertex_df[hit.vec,]
  
  # khroma::plot_scheme(khroma::color('BuRd')(7))
  max.logFC = max(abs(V(g)$logFC), na.rm = T)
  color_fun = circlize::colorRamp2(breaks = c(-max.logFC, -max.logFC/3*2, -max.logFC/3, 0,
                                              max.logFC/3, max.logFC/3*2, max.logFC),
                                   colors = khroma::color('BuRd')(7))
  V(g)$color = color_fun(V(g)$logFC)
  
  # V(g)$color = ifelse(!is.na(V(g)$hit.absmax), # if we have a hit
  #                     ifelse(V(g)$logFC>0, 'green', 'red'), # color green for up and red for down
  #                     'white') # otherwise keep white
  
  
  
  
  #### Plot the Graph ####
  
  pdf(file.path(out.dir,glue('{prefix}_network.pdf')), width = 16, height=16)
  plot.igraph(g, layout = matrix(c(V(g)$graphics_x, # coordinates matrix with x
                                   -V(g)$graphics_y), # flipped y (to match KEGG)
                                 ncol=2), # two columns
              main = title,
              vertex.label = V(g)$plot_name,
              # vertex.size = 15,
              vertex.size = ifelse(V(g)$type=="compound", 12, 15),
              # vertex.size = ifelse(V(g)$type=="compound", 3, 15*(nchar(V(g)$plot_name)/8)),
              # vertex.label.dist = ifelse(V(g)$type=="compound", .75, 0),
              vertex.size2 = 8,
              vertex.shape=V(g)$graphics_type)
  dev.off()
  
}



################################
####      Data Analysis     ####
################################

# helper function to determine if annot is blank / should be skipped
skip.annot = function(value_of_interest) {
  return( is.na(value_of_interest) ||
            value_of_interest=="" || # match blank
            grepl("^[nN]\\.*[aA]\\.*$",value_of_interest) || # match NA / N.A. of any case
            grepl("^[nN][aA][nN]$",value_of_interest)) # match NaN
}

# initialize log object to keep track of which annotations have what files
log_file = data.frame(annot.name = character(0),
                      subvalue.name = character(0), # valid unique subvalue name (for file lookup)
                      subvalue.value = character(0), # original subvaluename
                      valid.subvalue = logical(0),
                      valid.results = logical(0),
                      ## placeholder columns for significant features per test, and observed features across all pathways
                      N.Sign.DF = numeric(0),
                      N.Sign.DF.Compounds = numeric(0),
                      N.Sign.DF.Genes = numeric(0),
                      N.Obs.Overall = numeric(0),
                      N.Obs.Overall.Compounds = numeric(0),
                      N.Obs.Overall.Genes = numeric(0))
log_file_colsToDrop = c() # initialize empty vector to store columns to drop; to be populated in the loop below

cat("\n\n####################\nEnrichment Analysis\n\n")
for (annot_of_interest in names(annots)) {
  # skip annotation if it has too many values, or too few
  if (length(unique(annots[[annot_of_interest]])) > opt$max_annot_levels) { cat(glue("\n\n####################\nSkipping '{annot_of_interest}' annotation; too many annotation-values ({length(unique(annots[[annot_of_interest]]))}) to be considered discrete (>{opt$max_annot_levels}).\n\n")); next }
  if (length(unique(annots[[annot_of_interest]])) == 1 ) { cat(glue("\n\n####################\nSkipping '{annot_of_interest}' annotation; only one unique annotation-value.\n\n")); next }
  # otherwise run analysis for every unique annotation subvalue
  cat(glue("\n\n####################\nAnalyzing '{annot_of_interest}' Annotation \n####################\n\n"))
  
  #### Set Up Log-File / Directory ####
  # initialize log_file with all FALSE
  log_file.tmp = data.frame(annot.name = annot_of_interest,
                            subvalue.name = make.names(sort(unique(annots[[annot_of_interest]]))), # valid unique subvalue name (for file lookup)
                            subvalue.value = sort(unique(annots[[annot_of_interest]])), # original subvaluename
                            # use subvalue.name as row.names for quick lookup / edits
                            row.names = make.names(sort(unique(annots[[annot_of_interest]]))), # valid unique subvalue name (for file lookup)
                            valid.subvalue = FALSE,
                            valid.results = FALSE,
                            ## placeholder values for significant features per test, and observed features across all pathways
                            N.Sign.DF = NA_real_,
                            N.Sign.DF.Compounds = NA_real_,
                            N.Sign.DF.Genes = NA_real_,
                            N.Obs.Overall = NA_real_,
                            N.Obs.Overall.Compounds = NA_real_,
                            N.Obs.Overall.Genes = NA_real_)
  
  annot_dir = file.path(opt$output_dir, glue("results_{make.names(annot_of_interest)}"))
  dir.create(annot_dir)
  
  #### Run Enrichment on subvalues ####
  res.df.list = list() # initialize empty list for res.fin.df results
  network_annots[[annot_of_interest]] = list() # initialize empty list for pathway network-annotations
  for (value_of_interest.name in rownames(log_file.tmp)) {
    value_of_interest = log_file.tmp[value_of_interest.name,"subvalue.value"] # get original annot value for use in analysis
    if (skip.annot(value_of_interest)) { cat(glue("\nSkipping '{value_of_interest}' annotation.\n\n")); next }
    log_file.tmp[value_of_interest.name,"valid.subvalue"] = TRUE # mark subvalue as a valid subvalue
    
    ################################
    ####   Enrichment Analysis  ####
    ################################
    if(print_internal_placemarks) cat(glue("\n\n####################\nSetting Up {opt$anal_type} Analysis for {value_of_interest}\n####################\n\n"))
    #### Calculate Single-omic Enrichments ####
    if (opt$anal_type == "ORA") {
      if(print_internal_placemarks) cat("\n\n####################\nOverrepresenation Analysis on Metabolome\n\n")
      meta_ora = tryCatch(gct.to.ora.input(meta_input, annot_of_interest, value_of_interest, annots = annots,
                                           write_to_file = T, prefix = glue("{opt$output_prefix}_metabolome")),
                          error = function(e) {
                            if (e == "No residual degrees of freedom in linear model fits") { return(NULL) } else { stop(e) }
                          })
      if (is.null(meta_ora)) { cat(glue("\nSkipping '{value_of_interest}' annotation-value; linear model could not be fit for metabolome.\n\n")); next }
      meta_ora_vec = meta_ora$id
      res.meta = o.ea(meta_ora_vec, pathways,
                      background = background_meta,
                      min_overlap = opt$min_overlap)
      
      if (!multiomic && dim(res.meta$df)[1]==0) { cat(glue("\nSkipping '{value_of_interest}' annotation-value; no significant metabolite enrichments.\n\n")); next }
      
      # calculate genomic enrichment
      if (multiomic) {
        if(print_internal_placemarks) cat(glue("\n\n####################\nOverrepresenation Analysis on {opt$ome_type}\n\n"))
        ome_ora = tryCatch(gct.to.ora.input(ome_input, annot_of_interest, value_of_interest, annots = annots,
                                            write_to_file = T, prefix = glue("{opt$output_prefix}_{opt$ome_type}")),
                           error = function(e) {
                             if (e == "No residual degrees of freedom in linear model fits") { return(NULL) } else { stop(e) }
                           })
        if (is.null(ome_ora)) { cat(glue("\nSkipping '{value_of_interest}' annotation-value; linear model could not be fit for additional -ome.\n\n")); next }
        ome_ora_vec = ome_ora$id
        res.ome = o.ea(ome_ora_vec, pathways,
                       background = background_ome,
                       min_overlap = opt$min_overlap)
        if (dim(res.meta$df)[1]==0 && dim(res.ome$df)[1]==0) { cat(glue("\nSkipping '{value_of_interest}' annotation-value; no significant enrichments.\n\n")); next }
        
        if (opt$pval_comb=="query") {
          if(print_internal_placemarks) cat("\n\n####################\nOverrepresenation Analysis on Both Omes\n\n")
          res.mo = o.ea(c(meta_ora_vec, ome_ora_vec), pathways,
                        background = background_mo,
                        min_overlap = opt$min_overlap)
        }
      }
      # make logFC dataframe
      if (multiomic) {
        logFC_df = dplyr::select(rbind(meta_ora, ome_ora), c('id', 'logFC'))
      } else {
        logFC_df = dplyr::select(meta_ora, c('id', 'logFC'))
      }
      
      
    } else if (opt$anal_type == "QEA") {
      if(print_internal_placemarks) cat("\n\n####################\nQuantitative Enrichment Analysis on Metabolome\n\n")
      meta_qea = gct.to.qea.input(meta_input, annot_of_interest, value_of_interest, annots = annots,
                                  write_to_file = T, prefix=glue("{opt$output_prefix}_metabolome"))
      res.meta = q.ea(meta_qea$mat, meta_qea$cls, pathways, background_filter = opt$background_filter,
                      min_overlap = opt$min_overlap)
      if (!multiomic && dim(res.meta$df)[1]==0) { cat(glue("\nSkipping '{value_of_interest}' annotation-value; no significant metabolite enrichments.\n\n")); next }
      
      if (multiomic) {
        if(print_internal_placemarks) cat(glue("\n\n####################\nQuantitative Enrichment Analysis on {opt$ome_type}\n\n"))
        ome_qea = gct.to.qea.input(ome_input, annot_of_interest, value_of_interest, annots = annots,
                                   write_to_file = T, prefix=glue("{opt$output_prefix}_{opt$ome_type}"))
        res.ome = q.ea(ome_qea$mat, ome_qea$cls, pathways, background_filter = opt$background_filter,
                      min_overlap = opt$min_overlap)
        if (dim(res.meta$df)[1]==0 && dim(res.ome$df)[1]==0) { cat(glue("\nSkipping '{value_of_interest}' annotation-value; no significant enrichments.\n\n")); next }
        
        
        if (opt$pval_comb=="query") {
          if(print_internal_placemarks) cat("\n\n####################\nQuantitative Enrichment Analysis Both Omes\n\n")
          # wrangle data
          cid_intersect = intersect(meta_input@cid, ome_input@cid) # identify shared samples
          mo_input = merge_gct(meta_input, ome_input, dim = "row") %>% # merge GCT features
            subset_gct(cid = cid_intersect) # subset to shared samples
          mo_qea = gct.to.qea.input(mo_input, annot_of_interest, value_of_interest, annots = annots,
                                    write_to_file = T, prefix=glue("{opt$output_prefix}_multiomic"))
          # quantitative enrichment analysis
          res.mo = q.ea(mo_qea$mat, mo_qea$cls, pathways, background_filter = opt$background_filter,
                        min_overlap = opt$min_overlap)
        }
      }
      
    }
    
    #### Store phyper() parameters in log_file, drop irrelevant columns ####
    if (multiomic && opt$pval_comb=="query") {
      # combined analysis: store combined values in both columns
      if (opt$anal_type == "ORA") { # if ORA, store query size
        log_file.tmp[value_of_interest.name,"N.Sign.DF"] = res.mo$query.size
      } else { log_file_colsToDrop = c(log_file_colsToDrop, "N.Sign.DF") } # otherwise drop this column
      log_file.tmp[value_of_interest.name,"N.Obs.Overall"] = res.mo$universe.size
      log_file_colsToDrop = c(log_file_colsToDrop, "N.Sign.DF.Compounds", "N.Sign.DF.Genes", "N.Obs.Overall.Compounds", "N.Obs.Overall.Genes") # missing columns, to be dropped later
    } else if (multiomic) {
      # separate analyses: store separate values for compounds and genes
      if (opt$anal_type == "ORA") { # if ORA, store query size
        log_file.tmp[value_of_interest.name,"N.Sign.DF.Compounds"] = res.meta$query.size
        log_file.tmp[value_of_interest.name,"N.Sign.DF.Genes"] = res.ome$query.size
      } else { log_file_colsToDrop = c(log_file_colsToDrop, "N.Sign.DF.Compounds", "N.Sign.DF.Genes") } # otherwise drop these columns
      log_file.tmp[value_of_interest.name,"N.Obs.Overall.Compounds"] = res.meta$universe.size
      log_file.tmp[value_of_interest.name,"N.Obs.Overall.Genes"] = res.ome$universe.size
      log_file_colsToDrop = c(log_file_colsToDrop, "N.Sign.DF", "N.Obs.Overall") # missing columns, to be dropped later
    } else {
      # single-omic: store only compounds values, leave genes as NA
      if (opt$anal_type == "ORA") { # if ORA, store query size
        log_file.tmp[value_of_interest.name,"N.Sign.DF.Compounds"] = res.meta$query.size
      } else { log_file_colsToDrop = c(log_file_colsToDrop, "N.Sign.DF.Compounds") } # otherwise drop this column
      log_file.tmp[value_of_interest.name,"N.Obs.Overall.Compounds"] = res.meta$universe.size
      log_file_colsToDrop = c(log_file_colsToDrop, "N.Sign.DF", "N.Sign.DF.Genes", "N.Obs.Overall", "N.Obs.Overall.Genes") # missing columns, to be dropped later
    }
    
    log_file.tmp[value_of_interest.name,"valid.results"] = TRUE # mark subvalue as having enrichments
    
    #### Combine P-Values ####
    if (multiomic) {
      if (opt$pval_comb=="query") {
        res.df = res.mo$df
        res.fin.hits = res.mo$hits
      } else {
        res.meta.df = res.meta$df
        res.ome.df = res.ome$df
        all.paths = unique(union(rownames(res.meta.df),rownames(res.ome.df))) # get all paths with at least 1 hit in either ome
        shared.paths = intersect(rownames(res.meta.df),rownames(res.ome.df)) # get paths with hits in both omes
        
        if(print_internal_placemarks) cat("\n\n####################\nCombinging P-values\n\n")
        performWeightedZtest = function(p1, p2, w1, w2) {
          p.vec = c(p1,p2)
          weights = c(w1, w2)
          zp <- (qnorm(p.vec, lower.tail = FALSE) %*% weights)/sqrt(sum(weights^2))
          # res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE))
          pnorm(zp, lower.tail = FALSE)
        }
        if (opt$pval_comb=="pvalu") {
          w.m = 0.5
          w.g = 0.5
        } else if (opt$pval_comb=="pvalo") {
          # use universe sizes from result lists (filtered to dataset features)
          universe.meta = res.meta$universe.size # get universe size for metabolites
          universe.ome = res.ome$universe.size # get universe size for genes
          w.m <- universe.meta/(universe.meta+universe.ome)
          w.g <- universe.ome/(universe.meta+universe.ome)
        } else if (opt$pval_comb=="pvalp") {
          # use pathway entry counts from result dataframes (filtered to dataset features)
          pw.c = res.meta.df[all.paths,'Observed.Entries'] %>% ifelse(is.na(.), 0, .)
          pw.g = res.ome.df[all.paths,'Observed.Entries'] %>% ifelse(is.na(.), 0, .)
          w.m = pw.c/(pw.c+pw.g)
          w.g = pw.g/(pw.c+pw.g)
        } else stop(glue("Invalid p-value combination method '{opt$pval_comb}' selected. Please select either 'query' (combine queries), 'pvalu' (un-weighted), 'pvalo' (overall), or 'pvalp' (pathway-wise)."))
        
        # calculate combined p-values
        idx.to.comb = all.paths %in% shared.paths # identify pathways with hits in both datasets
        comb.pval = ifelse(idx.to.comb, # if the pathway had entries from both omes
                           # combine p-values
                           mapply(performWeightedZtest, # perform a weighted z-test to combine
                                  res.meta.df[all.paths,'Raw.P.Value'], # relevant p-value from compounds
                                  res.ome.df[all.paths,'Raw.P.Value'], # relevant p-value from genes
                                  w.m, w.g), # based on weights for compounds and genes respectively
                           # else select first non-missing
                           coalesce(as.vector(res.meta.df[all.paths,'Raw.P.Value']),
                                    as.vector(res.ome.df[all.paths,'Raw.P.Value'])))
        # create finalized data-frame, ordered by all.paths
        res.df = data.frame(ID = all.paths,
                            Pathway.Name = sapply(pathways, function(p) {p$name})[all.paths],
                            Test.Type = opt$anal_type) %>%
          { if (opt$anal_type == "ORA") {
            dplyr::mutate(.,
                          Observed.Entries = mapply(sum, res.meta.df[all.paths, "Observed.Entries"],
                                                         res.ome.df[all.paths, "Observed.Entries"],  # combined total
                                                    na.rm = TRUE), # ignoring NA values
                          N.Hits = mapply(sum, res.meta.df[all.paths, "N.Hits"],
                                               res.ome.df[all.paths, "N.Hits"], # combined total
                                          na.rm = TRUE),
                          # stats for genes
                          Observed.Entries.Genes = res.ome.df[all.paths, "Observed.Entries"],
                          N.Hits.Genes = res.ome.df[all.paths, "N.Hits"],
                          # stats for compounds
                          Observed.Entries.Compounds = res.meta.df[all.paths, "Observed.Entries"],
                          N.Hits.Compounds = res.meta.df[all.paths, "N.Hits"] )
          } else {
              dplyr::mutate(.,
                            Observed.Entries = mapply(sum, res.meta.df[all.paths, "Observed.Entries"],
                                                           res.ome.df[all.paths, "Observed.Entries"],  # combined total
                                                      na.rm = TRUE), # ignoring NA values
                            # don't create combined column for stats; this number would not be meaningful
                            # stats for genes
                            Observed.Entries.Genes = res.ome.df[all.paths, "Observed.Entries"],
                            Observed.Stat.Genes = res.ome.df[all.paths, "Observed.Stat"],
                            Expected.Stat.Genes = res.ome.df[all.paths, "Expected.Stat"],
                            # stats for compounds
                            Observed.Entries.Compounds = res.meta.df[all.paths, "Observed.Entries"],
                            Observed.Stat.Compounds = res.meta.df[all.paths, "Observed.Stat"],
                            Expected.Stat.Compounds = res.meta.df[all.paths, "Expected.Stat"])
          } } %>%
          # add combined p-value and adjust accordingly
          dplyr::mutate(Raw.P.Value = comb.pval,
                        neg.Log.P.Value = -log10(Raw.P.Value),
                        Holm.P.Value = p.adjust(Raw.P.Value, "holm"),
                        BH.P.Value = p.adjust(Raw.P.Value, "fdr")) %>%
          dplyr::arrange(Raw.P.Value) %>% # sort in order of significance
          dplyr::mutate_if(is.numeric, signif, 5) %>% # filter to 5 sigfigs
          column_to_rownames("ID") # add rownames back
        # merge hits into a single list, ordered by all.paths
        res.fin.hits = mapply(function(hits.m, hits.g) {c(as.vector(hits.m), as.vector(hits.g))}, # force vector format to avoid empty lists if one dataset is empty
                              res.meta$hits[all.paths], res.ome$hits[all.paths], SIMPLIFY = F)
      }
    } else {
      res.df = res.meta$df %>%
        rownames_to_column("ID") %>% # save rownames in a column
        dplyr::mutate_if(is.numeric, signif, 5) %>% # filter to 5 sigfigs
        column_to_rownames("ID") # add rownames back
      res.fin.hits = res.meta$hits
    }
    
    
    #### Wrangle Pathway-Network Annotations into a list ####
    if(print_internal_placemarks) cat("\n\n####################\nWrangling Pathway Network-Annotations\n\n")
    network_annots_tmp = list() # initialize empty list for this annotation's network-annotations
    if (exists("pthw_json")) { # if we have pathway information
      for (pthw in names(pthw_json)) {
        network_annots_pthw_tmp = dplyr::select(pthw_json[[pthw]]$nodes, c(id, names))
        network_annots_pthw_tmp$entries = sapply(network_annots_pthw_tmp$names, function (names) {
          strsplit(names, " ", fixed = TRUE) %>% unlist() # unlist entries
        })
        network_annots_pthw_tmp$hits.all = sapply(network_annots_pthw_tmp$entries, function (entries) {
          hits.sign = all_features[which(all_features %in% entries)] %>% as.character() 
          return(hits.sign)
        })
        network_annots_pthw_tmp$hits.all.hr = sapply(network_annots_pthw_tmp$hits.all, function (ft) {
          unname(feature_map[ft]) %>% as.character() 
        })
        ## if we have logFC data, get the significant features + corresponding logFC values
        if(exists("logFC_df")) {
          # add all significant hits
          network_annots_pthw_tmp$hits.sign = sapply(network_annots_pthw_tmp$entries, function (entries) {
            logFC_df$id[which(logFC_df$id %in% entries)] %>% as.character() 
          })
          network_annots_pthw_tmp$hits.sign.hr = sapply(network_annots_pthw_tmp$hits.sign, function (ft) {
            unname(feature_map[ft]) %>% as.character() 
          })
          # add logFC of all significant hits
          network_annots_pthw_tmp$logFC.all.sign = sapply(network_annots_pthw_tmp$entries, function (entries) {
            logFC_df$logFC[which(logFC_df$id %in% entries)] %>% as.numeric() 
          })
          # add feature ID with highest logFC
          network_annots_pthw_tmp$hit.absmax = sapply(network_annots_pthw_tmp$entries, function (entries) {
            logFC.tmp = logFC_df[which(logFC_df$id %in% entries),]
            logFC = logFC.tmp[which.max(abs(logFC.tmp$logFC)), 'id']
            return(as.character(logFC))
          })
          network_annots_pthw_tmp$hit.absmax.hr = sapply(network_annots_pthw_tmp$hit.absmax, function (ft) {
            unname(feature_map[ft]) %>% as.character() 
          })
          # add abs.max logFC
          network_annots_pthw_tmp$logFC.absmax = sapply(network_annots_pthw_tmp$entries, function (entries) {
            logFC.tmp = logFC_df[which(logFC_df$id %in% entries),]
            logFC = logFC.tmp[which.max(abs(logFC.tmp$logFC)), 'logFC']
            if(length(logFC)==0) logFC = NA # add NA placeholder if there were no logFC values
            return(as.numeric(logFC))
          })
        } else {
          network_annots_pthw_tmp = network_annots_pthw_tmp %>%
            mutate(hits.sign = "no ORA performed",
                   hits.sign.hr = NA,
                   logFC.all.sign = NA,
                   hit.absmax = NA,
                   hit.absmax.hr = NA,
                   logFC.absmax = NA)
        }
        
        network_annots_tmp[[pthw]] = network_annots_pthw_tmp
        
      }
    }
    network_annots[[annot_of_interest]][[value_of_interest.name]] = network_annots_tmp
    
    #### Topological Measurements ####
    if (has_topology) {
      if(print_internal_placemarks) cat("\n\n####################\nCalculating Topological Scores\n\n")
      calculate.impact = function(pathway.hits, pathway_topology_scores,
                                  score_type = "bc") {
        if (! glue("{score_type}.list") %in% names(pathway_topology_scores[[1]])) stop(glue("Invalid score-type {score_type} chosen. Please select from 'bc', 'cc', or 'dc', and ensure that pathways contain topological scores metrics."))
        
        impact.vec = sapply(names(pathway.hits), function(p) { # for each pathway
          # boolean vector of grouped-entries with hits
          entries.sep = sapply(pathway_topology_scores[[p]]$entries.grouped, function(x) {strsplit(x, " ", fixed = TRUE)} )
          hit.vec = sapply(entries.sep, function(entr) {
            any(pathway.hits[[p]] %in% entr)
          })
          # take sum of topological measure for the relevant pathway
          sum(pathway_topology_scores[[p]][[glue("{score_type}.list")]][hit.vec])
        })
        return(impact.vec)
      }
      
      impact.bc <- calculate.impact(res.fin.hits, pathway_topology_scores, 'bc')
      impact.cc <- calculate.impact(res.fin.hits, pathway_topology_scores, 'cc')
      impact.dc <- calculate.impact(res.fin.hits, pathway_topology_scores, 'dc')
      
      res.fin.df = res.df %>%
        rownames_to_column('ID') %>%
        dplyr::mutate(Impact.BC = impact.bc[ID],
                      Impact.CC = impact.cc[ID],
                      Impact.DC = impact.dc[ID]) %>%
        dplyr::mutate_if(is.numeric, signif, 5) %>% # filter to 5 sigfigs
        column_to_rownames('ID')
    } else { res.fin.df = res.df }
    
    
    #### Exporting Results ####
    res.df.list[[value_of_interest.name]] = res.fin.df
    
    fn = glue("{opt$output_prefix}_{opt$anal_type}_metabolome_")
    if (multiomic) fn = paste0(fn, glue("with.{opt$ome_type}_by.{opt$pval_comb}_")) # add multiomic info
    fn = paste0(fn, glue("{make.names(annot_of_interest)}_{value_of_interest.name}_results.csv"))
    write.csv(res.fin.df, file = file.path(annot_dir,fn))
    
    
    
    ################################
    ####   Figure Generation    ####
    ################################
    if(print_internal_placemarks) cat("\n\n####################\nFigure Generation\n\n")
    
    #### logP vs Impact ####
    if (has_topology) {
      # make 'volcano' plot
      ggplot(res.fin.df, aes(x = !!opt$impact_metric, y = neg.Log.P.Value, color = neg.Log.P.Value)) + 
        geom_point() +
        # xlab("Impact (Degree Centrality)") +
        ylab("-log(P-Value)")
      # save to file
      fn = glue("{opt$output_prefix}_{make.names(annot_of_interest)}_{value_of_interest.name}_negLogFC_vs_{opt$impact_metric}.png")
      ggsave(file.path(annot_dir, fn))
    }
    
    #### Enrichment Ratio-- % Pathway Hits vs Significance ####
    if (opt$anal_type == "ORA") { # only plot enrichment ratio for ORA
      plot_df = head(res.fin.df,opt$top_n_networks) %>%
        dplyr::arrange(neg.Log.P.Value)
      # make 'volcano' plot
      ggplot(plot_df, aes(x = neg.Log.P.Value,
                          y = factor(Pathway.Name, level = Pathway.Name))) + 
        geom_point(aes(size = Observed.Entries/Observed.Entries), color = 'black', shape=1) + # create outline showing 100%
        geom_point(aes(size = N.Hits/Observed.Entries), color = 'red') + # create inner dot showing the % of hits
        scale_size_continuous(labels = scales::percent,
                              name = 'Pathway Coverage\n(# Hits / Pathway Entries)') +
        ylab("Pathway") +
        xlab("-log(P.Value)")
        
      # save to file
      fn = glue("{opt$output_prefix}_{make.names(annot_of_interest)}_{value_of_interest.name}_EnrichemntRatio.png")
      ggsave(file.path(annot_dir, fn))
      # file.copy(file.path(annot_dir, fn), '/opt/input/', overwrite=T)
    }

    #### Network Graphs ####
    if (!is.null(graph_list) & exists("logFC_df")) {
      val_dir = file.path(annot_dir, glue("networks_{value_of_interest.name}"))
      dir.create(val_dir)
      
      for (p in rownames(head(res.fin.df,opt$top_n_networks))) { # create network plot for top opt$top_n_networks networks
        plot.network(g = graph_list[[p]], hits = res.fin.hits[[p]], logFC_df = logFC_df,
                     title = res.fin.df[p, 'Pathway.Name'],
                     write_to_file = T, prefix=glue("{opt$output_prefix}_{make.names(annot_of_interest)}_{value_of_interest.name}_{p}"),
                     out.dir = val_dir)
      }
    }
  }
  
  #### Add to Log File (even if there was no enrichment) ####
  rownames(log_file.tmp) = NULL # drop rownames
  log_file = rbind(log_file, log_file.tmp) # append tmp logfile to logfile
  
  ################################
  #### Compile Results ####
  ################################
  # skip compilation if nothing had enrichment
  if (length(res.df.list)==0) { cat(glue("\n\n####################\nNo significant enrichments found for '{annot_of_interest}' annotation.\n\n")); next }
  
  #### Exporting Full-Annot Excel ####
  fn = glue("{opt$output_prefix}_{opt$anal_type}_metabolome_")
  if (multiomic) fn = paste0(fn, glue("with.{opt$ome_type}_by.{opt$pval_comb}_")) # add multiomic info
  fn = paste0(fn, glue("{make.names(annot_of_interest)}_fullResults.xls"))
  
  WriteXLS::WriteXLS(res.df.list,
                     ExcelFileName=file.path(opt$output_dir,fn),
                     FreezeRow=1, FreezeCol=1,
                     SheetNames=make.unique(substr(names(res.df.list), 1, 20)),
                     row.names=T, BoldHeaderRow=T, AutoFilter=T)
  
  
  #### Heatmap of -log(P.Value) for Top Pathways ####
  res.df.long = lapply(1:length(res.df.list), function(i) {
    res.df = res.df.list[[i]]
    df = dplyr::mutate(res.df, Annot.Subvalue = names(res.df.list)[i])
  }) %>% do.call(rbind,.)
  
  sign_pathways = filter(res.df.long, BH.P.Value<opt$pval_signif)$Pathway.Name
  if (length(sign_pathways) > max_hm_pathways) {
    top_n_pathways = max(1,floor(max_hm_pathways/length(res.df.list))) # pathways to plot per subvalue (minimum 1 per subvalue)
    sign_pathways = lapply(res.df.list, function(res.df) {
      df = dplyr::filter(res.df, BH.P.Value<opt$pval_signif) %>% # filter to significant values
        slice_min(order_by = BH.P.Value, # order by corrected p-value
                  n = top_n_pathways, # select n features with lowest p-values
                  with_ties = TRUE) # keep ties
      return(df$Pathway.Name)
    }) %>% { unique(unlist(.)) } # get only the unique elements
  }
  
  value_col = 'neg.Log.P.Value'
  heatmap_df = dplyr::select(res.df.long, c('Pathway.Name', !!value_col, 'Annot.Subvalue')) %>%
    pivot_wider(id_cols = c('Pathway.Name'), names_from = 'Annot.Subvalue',
                values_from = value_col) %>%
    filter(Pathway.Name %in% sign_pathways) %>%
    column_to_rownames('Pathway.Name') %>%
    replace(is.na(.), 0) # replace NA -log(P.Value) with 0 for clustering
  
  sign_df = dplyr::select(res.df.long, c('Pathway.Name', 'BH.P.Value', 'Annot.Subvalue')) %>%
    pivot_wider(id_cols = c('Pathway.Name'), names_from = 'Annot.Subvalue',
                values_from = 'BH.P.Value') %>%
    filter(Pathway.Name %in% sign_pathways) %>%
    column_to_rownames('Pathway.Name')
  
  # generate heatmap
  hm.title = glue("{opt$output_prefix} {opt$anal_type} Results\nTop Significant Pathways for {annot_of_interest}")
  # hm.subtitle = glue("{value_col}")
  hm <- Heatmap(as.matrix(heatmap_df), # plot normalized heatmap (as matrix to suppress warning)
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(is.na(sign_df[i, j])) {
                    grid.text("", x, y)
                  # } else if(sign_df[i, j] <0.01) {
                  #   grid.text("***", x, y, gp = gpar(fontface = "bold", col = "white"))
                  # } else if(sign_df[i, j] <0.02) {
                  #   grid.text("**", x, y, gp = gpar(fontface = "bold", col = "white"))
                  } else if(sign_df[i, j] <opt$pval_signif ) {
                    grid.text("*", x, y, gp = gpar(fontface = "bold", col = "white"))
                  }
                }, # add significance stars
                col = circlize::colorRamp2(breaks = c(0,seq(1e-10, 5, length.out=7)),
                                           colors = c('grey',khroma::color('YlOrBr')(7))), # grey for 0 (NA), YlOrBr range for 0->5
                width = unit(dim(heatmap_df)[2]*8, "mm"),
                height = unit(dim(heatmap_df)[1]*8, "mm"),
                cluster_rows = TRUE, row_dend_side = "right",
                cluster_columns = FALSE, column_dend_side = "bottom",
                show_row_names = TRUE, row_names_side = "left",
                column_title = hm.title,
                column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                # column_title = hm.subtitle,
                # column_title_gp = gpar(fontsize = 16),
                show_heatmap_legend = T,
                row_title_rot = 0, # horizontal titles
                heatmap_legend_param = list(title = value_col, legend_width = unit(30, "mm"),
                                            direction = "horizontal"))

  # draw heatmap to PDF
  fn = glue("{opt$output_prefix}_{opt$anal_type}_metabolome_")
  if (multiomic) fn = paste0(fn, glue("with.{opt$ome_type}_by.{opt$pval_comb}_")) # add multiomic info
  fn = paste0(fn, glue("{make.names(annot_of_interest)}_heatmap"))
  pdf(file.path(opt$output_dir,paste0(fn,'.pdf')),
      width = dim(heatmap_df)[2]*8/25.4+8, # adapt to width of heatmap + 4
      height = dim(heatmap_df)[1]*8/25.4+4) # adapt to height of heatmap + 2 inches for header/footer
  draw(hm, padding = unit(c(4, 4, 4, 4), "mm"),
       # column_title = hm.title, # add title
       # column_title_gp = gpar(fontsize = 16, fontface = "bold"), #bold and increase font
       heatmap_legend_side = 'bottom')
  dev.off()
  # and export to PNG
  png(file.path(opt$output_dir,paste0(fn,'.png')), units="in", res=300,
      width = dim(heatmap_df)[2]*8/25.4+8, # adapt to width of heatmap + 4
      height = dim(heatmap_df)[1]*8/25.4+4) # adapt to height of heatmap + 2 inches for header/footer
  draw(hm, padding = unit(c(4, 4, 4, 4), "mm"),
       # column_title = hm.title, # add title
       # column_title_gp = gpar(fontsize = 16, fontface = "bold"), #bold and increase font
       heatmap_legend_side = 'bottom')
  dev.off()
  
  #### more figures...? ####
  
}

# write log file
log_file = dplyr::select(log_file, -log_file_colsToDrop ) # drop all columns that only contain NA values
write.csv(log_file,
          file = file.path(opt$output_dir, glue("{opt$output_prefix}_log_file.csv")))
# write pathway network-annotations to JSON
exportJSON <- toJSON(network_annots)
write_json(exportJSON, path=file.path(opt$output_dir, glue("{opt$output_prefix}_network_annots.JSON")))

################################
####        Tar File        ####
################################
cat(glue("\n\n####################\nTarring Files\n\n"))
fn = paste0(opt$output_prefix,'_MetaboAnalyst.tar.gz')

pwd = getwd()
setwd(opt$output_dir)
tar(tarfile = fn, compression = "gzip", tar="tar") # compress files in output directory
file.copy(fn, pwd, overwrite = T)
setwd(pwd)





