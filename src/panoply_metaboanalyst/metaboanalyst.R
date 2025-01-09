#!/usr/bin/env Rscript
#
# Copyright (c) 2024 The Broad Institute, Inc. All rights reserved.
#
rm(list=ls())
options( warn = -1, stringsAsFactors = F )
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
  # make_option( c("-w", "--pathway_gmt"), action='store', type='character',  dest='pathway_gmt', help='GMT file containing pathways of interest.'),
  make_option( c("-g", "--groups_file"), action='store', type='character',  dest='groups_file', help='Groups-file, i.e. an annotations file subsetted to annotations of interest.'),
  #### Analysis ####
  make_option( c("-l", "--max_annot_levels"), action='store', type='numeric', dest='max_annot_levels', help='Maximum number of levels an annotation can have and be considered discrete.'), # default='10'),
  make_option( c("-a", "--anal_type"), action='store', type='character', dest='anal_type', help='Analysis method to use ("ORA" for Overrepresentation Analysis or "QEA" for Quantitative Enrichment Analysis).'), 
  make_option( c("-b", "--pval_comb"), action='store', type='character', dest='pval_comb', help='Method for combining p-values in multiomic enrichment analysis. Options include "query" (combine queries), "pvalu" (unweighted), "pvalo" (overall), or "pvalp" (pathway-level).'), 
  make_option( c("-p", "--pval_signif"), action='store', type='numeric', dest='pval_signif', help='P-value threshold for significant enrichement.'), 
  make_option( c("-k", "--top_n_networks"), action='store', type='numeric', dest='top_n_networks', help='Top N networks to plot per annot subvalue.'), 
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.'),  # default = 2),
  make_option( c("-f", "--output_directory"), action='store', type='character',  dest='output_dir', help='Directory to output files to.',  default = 'results/'),
  make_option( c("-y", "--yaml"), action="store", dest='yaml_file', type="character", help="Path to .yaml file with parameters."),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)

#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                       # # for testing arguments
                       # args = c(
                       #   # '--metabolome_gct',"/opt/input/HMDB_ID_GCTs/ODG-v2_2-metabolomics_log_norm-HMDB_UNIQUE.gct",
                       #   '--metabolome_gct',"/opt/input/ODG-v3-metabolome-polar-log2-median-norm-QCfilter.gct",
                       #   '-n',"hmdb_id",
                       #   '-i',"HMDB.ID",
                       #   # '-n',"kegg_id",
                       #   # '-i',"KEGG.ID",
                       #   '--ome_gct',"/opt/input/ODG-v3-proteome-SpectrumMill-ratio-QCfilter-NArm.gct",
                       #   '-t',"prot",
                       # #   # '--ome_gct',"opt/input/ODG-v2_2-rnaseq-expression-TPM-protein-coding-log2-median-norm-NArm-with-NMF.gct",
                       # #   # '-t',"RNA",
                       # #   # '-g',"opt/input/sample-info.csv",
                       #   '-g',"opt/input/groups-subset.csv",
                       # #   '-l',"15",
                       # #   '-a',"QEA",
                       # #   '-b',"pvalp",
                       # #   '-p',"0.01",
                       # #   '-k',"15",
                       #   '-y',"opt/input/master-parameters.yaml",
                       # #   # '-f',"/opt/input/prelim_results",
                       #   '-x',"ODG_v3")
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

compound_map = qs::qread(file.path(opt$lib_dir, "pathway_db/compound_db.qs")) # mapping between ID types for metabolic compounds


#### Read in / Format Pathways ####
if (multiomic) {
  # format QS object from MetaboAnalyst Site into GMT-adjacent
  pathways_qs = qs::qread(file.path(opt$lib_dir,"pathway_db/hsa.qs"))
  pathways = lapply(unname(pathways_qs$path.ids), function(pathway) {
    pathway_list = list(ID = pathway,
                        name = names(pathways_qs$path.ids[pathways_qs$path.ids==pathway]),
                        entries = unlist(sapply(pathways_qs$mset.list[[pathway]], USE.NAMES=F, # unlist pathway elements
                                                function(x) {strsplit(x, " ", fixed = TRUE)} )), # and separate strings into single entries
                        cmpd.counts = sum(grepl(glue("^(cpd)|(gl):"), pathways_qs$mset.list[[pathway]])), # count compounds (before unlisting)
                        gene.counts = sum(grepl(glue("^hsa:"), pathways_qs$mset.list[[pathway]]))) # count genes (before unlisting)
    return(pathway_list)
  })
  names(pathways) = pathways_qs$path.ids
  # set pathway IDs
  pathway_id_type = list(meta="kegg_id",
                         gene="ENTREZID")
  
  # get unique compound / gene counts
  unique.cmpd = pathways_qs$uniq.cmpd.count
  unique.gene = pathways_qs$uniq.gene.count
  
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
  pathways_gmt = read.GMT(file.path(opt$lib_dir, "pathway_db/smpdb_pathway.gmt"))
  pathways = lapply(pathways_gmt, function(p) {
    pathway_list = list(ID = p$id,
                        name = p$name,
                        entries = p$genes,
                        cmpd.counts = length(p$genes))
    return(pathway_list)
  })
  # set pathway IDs
  pathway_id_type = list(meta="hmdb_id")
  
  # get unique compound counts
  unique.cmpd = length(unique(unlist(lapply(pathways_gmt, function(p) {p$genes})))) # get unique genes
  
  # binary toggle for topological measures
  has_topology = FALSE
  # graph list
  graph_list = NULL
}


#### Write Networks to JSON ####
if (!is.null(graph_list)) {
  pthw_json = list()
  for (pthw in names(graph_list)) {
    pthw_json[[pthw]] = d3r::d3_igraph(igrf = graph_list[[pthw]], json = FALSE)
  }
  
  exportJSON <- toJSON(pthw_json)
  write_json(exportJSON, path="networks.JSON")
  file.copy('networks.JSON', opt$output_dir) # copy network JSON to output directory
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

valid_cpd_rid = gct_meta@rid[which(cpd_vec %in% compound_map[[opt$meta_id_type]])]
if (length(valid_cpd_rid)==0) stop(glue("No IDs in the GCT mapped to valid compounds. Please check that your data uses {opt$meta_id_type} IDs, or select a different ID type."))
# subset to valid compound IDs
# toDo: add lipid ID mapping
cat(glue("\nOut of {length(cpd_vec)} features, {length(valid_cpd_rid)} mapped to valid compound IDs.\n\n"))
gct_meta_filt = subset_gct(gct_meta, rid=valid_cpd_rid)
# map to new ID type if necessary
if (opt$meta_id_type != pathway_id_type$meta || # if we need to change ID type
    !is.null(opt$meta_id_col)) { # OR the id column wasn't the rid
  # get relevant IDs from rid or rdesc
  if (is.null(opt$meta_id_col)) { cpd_vec_filt = gct_meta_filt@rid } else { cpd_vec_filt = gct_meta_filt@rdesc[[opt$meta_id_col]] }
  new_rid = compound_map[match(cpd_vec_filt, compound_map[[opt$meta_id_type]]), # match subsetted GCT rid to compound map
                         pathway_id_type$meta] # overwrite with new ID type
  rid_dup = !(new_rid %in% unique(new_rid[duplicated(new_rid)])) # identify duplicated RIDs (covers NA values)
  # overwrite RID and drop duplicated
  tmp = overwrite_rid(gct_meta_filt, new_rid, allow_dups = T) 
  meta_val = subset_gct(tmp, rid_dup)
  cat(glue("\nOut of {length(cpd_vec_filt)} features, {length(meta_val@rid)} had valid and unique {pathway_id_type$meta} IDs.\n\n"))
  
  # map between IDs and human-readable names
  feature_map = compound_map$name # get compound names
  names(feature_map) = compound_map[, c(pathway_id_type$meta)] # name with pathway IDs
} else { meta_val = gct_meta_filt }



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
  meta_mo = overwrite_rid(meta_val, paste0("cpd:", meta_val@rid))
  ome_mo = overwrite_rid(ome_val, paste0("hsa:", ome_val@rid))
  
  # map between IDs and human-readable names
  gene_sym_map_df = AnnotationDbi::select(org.Hs.eg.db, ome_val@rid, 'SYMBOL', pathway_id_type$gene) # get map between pathway IDs and geneSymbols
  gene_sym_map = gene_sym_map_df$SYMBOL # get geneSymbols
  names(gene_sym_map) = paste0("hsa:", gene_sym_map_df[,pathway_id_type$gene]) # match with pathway IDs (with hsa: prefix)
  # append gene IDs to feature_map
  names(feature_map) = paste0("cpd:", names(feature_map)) # add cpd prefix to feature_map
  feature_map = c(feature_map, gene_sym_map)
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


################################
####  Perform Enrichement   ####
################################

#### Quantitative Enrichment Analysis Function ####
# takes GCT & annot of interest
# returns transposed matrix & binary numeric cls vector
gct.to.qea.input = function(gct, annot_of_interest, value_of_interest, annots = NULL,
                                 write_to_file = F, prefix="results") {
  #### Data Wrangling ####
  if (is.null(annots)) { annots_fin = gct@cdesc } else {
    if (sum(gct@cid %in% annots[['Sample.ID']]) == 0 ) stop("No samples had matching IDs in provided annotation table.")  # check that we have overlapping IDs
    annots_fin = annots[match(gct@cid, annots[['Sample.ID']]),] # reorder to match gct@cid
  }
  # create binary CLS
  cls = annots_fin[[annot_of_interest]] %>%
    { ifelse(!is.na(.) & .==value_of_interest, 1,0) } # NA -> no_<value_of_interest>
  # toDo: consider whether we want to include NA as "not", or drop them entirely
  
  # pull out matrix and transpose so columns are features
  mat = as.matrix(t(gct@mat))
  # optionally save input dataset to a file
  if (write_to_file) {
    write.csv(data.frame(cls = cls, mat),
              file=glue("{prefix}_QEA_input_{make.names(annot_of_interest)}_{value_of_interest.name}_vs_{value_of_interest.name}_not.csv"))
  }
  return(list(mat = mat,
              cls = cls))
}
q.ea = function(mat, cls, pathways, uniq.len = NULL, p.val.min=2.3233E-11) {
  if (dim(mat)[1] != length(cls)) stop(glue("Matrix ({dim(mat)[1]}) and CLS vector ({length(cls)}) have differing dimensions. Please ensure that your matrix has samples for rows and columns for features, and that your CLS has an entry for every sample."))
  if (!is.numeric(cls)) stop(glue("CLS vector is not numeric! Please supply a binary numeric vector to annotate your samples."))
  
  #### Enrichment Analysis ####
  # number of entries in pathways
  if(is.null(uniq.len)) uniq.len = sapply(pathways, function(p) {sum(length(unique(p$entries)))})
  # number of hits in pathways
  pathway.hits = sapply(pathways, function(p) {
    x = p$entries
    unique(x[x %in% colnames(mat)])
  })
  if ( max(sapply(pathway.hits, length))==0 ) stop("No hits found in any pathway.")
  
  # calculate enrichment
  gt.obj <- globaltest::gt(cls, mat, subsets = pathway.hits)
  gt.res <- globaltest::result(gt.obj)
  
  # wrangle outputs
  res.df = data.frame(ID = names(pathways),
                      Pathway.Name = unname(sapply(pathways, USE.NAMES=F, function(p) {p$name})), # pull name into DF
                      N.Entries = uniq.len,
                      N.Hits = sapply(pathway.hits, length)) %>%
    dplyr::mutate(Observed = gt.res[, 2],
                  Expected = gt.res[, 3],
                  Enrichment.Ratio = Observed/Expected) %>% # Enrichment.Ratio
    dplyr::mutate(Test.Type = "QEA",
                  Raw.P.Value =  gt.res[, 1] %>% ifelse(.==0, p.val.min, .), # add a minimum p-value of 2.3233E-11, to avoid -log(P.Value)
                  neg.Log.P.Value = -log10(Raw.P.Value),
                  Holm.P.Value = p.adjust(Raw.P.Value, "holm"),
                  BH.P.Value = p.adjust(Raw.P.Value, "fdr")) %>%
    dplyr::filter(N.Hits>0) %>% # filter to valid hits
    dplyr::arrange(Raw.P.Value) %>% # sort in order of significance
    # dplyr::mutate_if(is.numeric, signif, 5) %>% # filter to 5 sigfigs # DONT round before combination
    column_to_rownames("ID") # add rownames back
  
  return(list(df = res.df,
              hits = pathway.hits))                  
}


#### Hypergeometric Function ####
# takes GCT & annot of interest
# returns DF with significant features and LogFC values
gct.to.ora.input = function(gct, annot_of_interest, value_of_interest, annots = NULL,
                            significance_cutoff = 0.05,
                            write_to_file = F, prefix="results") {
  #### Data Wrangling ####
  if (is.null(annots)) { annots_fin = gct@cdesc } else {
    if (sum(gct@cid %in% annots[['Sample.ID']]) == 0 ) stop("No samples had matching IDs in provided annotation table.")  # check that we have overlapping IDs
    annots_fin = annots[match(gct@cid, annots[['Sample.ID']]),] # reorder to match gct@cid
  }
  # create binary CLS
  cls = annots_fin[[annot_of_interest]] %>%
    { ifelse(!is.na(.) & .==value_of_interest, ., glue("{value_of_interest}_not")) } # NA -> no_<value_of_interest>
  # toDo: consider whether we want to include NA as "not", or drop them entirely
  
  #### T-Test ####
  require(limma)
  source('https://raw.githubusercontent.com/broadinstitute/protigy/master/src/modT.R') # for modT.test.2class()
  d = rownames_to_column(as.data.frame(gct@mat), "feature_id")
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
o.ea = function(queries, pathways, uniq.len=NULL, uniq.count=NULL, p.val.min=2.3233E-11) {
  # calculate analysis metrics to be used later
  current.universe <- unique(unlist(lapply(pathways, function(p) {p$entries}))) # get all unique entries in pathway
  # if we aren't limiting the feature-space (i.e. only hsa: or only cpd: / gl:), calculate size of featurespace
  if(is.null(uniq.len)) uniq.len = sapply(pathways, function(p) {sum(length(unique(p$entries)))})
  if(is.null(uniq.count)) uniq.count = sum(length(current.universe))
  # calculate hits
  queries.subset <- queries[queries %in% current.universe] # subset to queries in current.universe
  pathway.hits <- lapply(pathways, function(p) { # for each pathway
    y = queries.subset %in% unlist(p$entries) # identify which queries are hits
    queries.subset[y] # return those hits
  })
  pathway.hits.num = sapply(pathway.hits, length) # count the number of hits per pathway
  
  # perform hypergeometric test (equal end of fisher exact test)
  p.val = phyper(pathway.hits.num - 1, uniq.len, uniq.count - uniq.len,
                 length(queries.subset), lower.tail = F) %>%
    ifelse(.==0, p.val.min, .) # add a minimum p-value (2.3233E-11), to avoid -log(P.Value)=Inf
  
  res.df = data.frame(ID = sapply(pathways, function(p) {p$ID}),
                      Pathway.Name = sapply(pathways, function(p) {p$name}),
                      N.Entries = uniq.len,
                      N.Hits = pathway.hits.num) %>%
    dplyr::mutate(Observed = N.Hits,
                  Expected = length(queries.subset) * (uniq.len/uniq.count),
                  Enrichment.Ratio = Observed/Expected) %>% # Enrichment.Ratio
    dplyr::mutate(Test.Type = "ORA",
                  Raw.P.Value = p.val,
                  neg.Log.P.Value = -log10(Raw.P.Value),
                  Holm.P.Value = p.adjust(Raw.P.Value, "holm"),
                  BH.P.Value = p.adjust(Raw.P.Value, "fdr")) %>%
    dplyr::filter(N.Hits>0) %>% # filter to valid hits
    dplyr::arrange(Raw.P.Value) %>% # sort in order of significance
    # dplyr::mutate_if(is.numeric, signif, 5) %>% # filter to 5 sigfigs # DONT round before combination
    column_to_rownames('ID')
  
  return(list(df = res.df,
              hits = pathway.hits))
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

# initialize log object to keep track of which annotations have what files
log_file = data.frame(annot.name = character(0),
                      subvalue.name = character(0), # valid unique subvalue name (for file lookup)
                      subvalue.value = character(0), # original subvaluename
                      valid.subvalue = logical(0),
                      valid.results = logical(0))

cat("\n\n####################\nEnrichment Analysis\n\n")
for (annot_of_interest in names(annots)) {
  # skip annotation if it has too many values, or too few
  if (length(unique(annots[[annot_of_interest]])) > opt$max_annot_levels) { cat(glue("\n\n####################\nSkipping '{annot_of_interest}' annotation; too many annotation-values ({length(unique(annots[[annot_of_interest]]))}) to be considered discrete (>{opt$max_annot_levels}).\n\n")); next }
  if (length(unique(annots[[annot_of_interest]]))  == 1 ) { cat(glue("\n\n####################\nSkipping '{annot_of_interest}' annotation; only one unique annotation-value.\n\n")); next }
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
                            valid.results = FALSE)
  
  annot_dir = file.path(opt$output_dir, glue("results_{make.names(annot_of_interest)}"))
  dir.create(annot_dir)
  
  #### Run Enrichment on subvalues ####
  res.df.list = list() # initialize empty list for res.fin.df results
  network_annots[[annot_of_interest]] = list() # initialize empty list for pathway network-annotations
  for (value_of_interest.name in rownames(log_file.tmp)) {
    value_of_interest = log_file.tmp[value_of_interest.name,"subvalue.value"] # get original annot value for use in analysis
    if (is.na(value_of_interest) || value_of_interest=="") { cat(glue("\nSkipping '{value_of_interest}' annotation.\n\n")); next }
    log_file.tmp[value_of_interest.name,"valid.subvalue"] = TRUE # mark subvalue as a valid subvalue
    
    ################################
    ####   Enrichment Analysis  ####
    ################################
    if(print_internal_placemarks) cat(glue("\n\n####################\nSetting Up {opt$anal_type} Analysis for {value_of_interest}\n####################\n\n"))
    #### Calculate Single-omic Enrichments ####
    if (opt$anal_type == "ORA") {
      if(print_internal_placemarks) cat("\n\n####################\nOverrepresenation Analysis on Metabolome\n\n")
      meta_ora = gct.to.ora.input(meta_input, annot_of_interest, value_of_interest, annots = annots,
                                  write_to_file = T, prefix = glue("{opt$output_prefix}_metabolome"))
      meta_ora_vec = meta_ora$id
      res.meta = o.ea(meta_ora_vec, pathways,
                      uniq.count = unique.cmpd,
                      uniq.len = sapply(pathways, function(p) {p$cmpd.counts}))
      
      if (!multiomic && dim(res.meta$df)[1]==0) { cat(glue("\nSkipping '{value_of_interest}' annotation-value; no significant metabolite enrichments.\n\n")); next }
      
      # calculate genomic enrichment
      if (multiomic) {
        if(print_internal_placemarks) cat(glue("\n\n####################\nOverrepresenation Analysis on {opt$ome_type}\n\n"))
        ome_ora = gct.to.ora.input(ome_input, annot_of_interest, value_of_interest, annots = annots,
                                   write_to_file = T, prefix = glue("{opt$output_prefix}_{opt$ome_type}"))
        ome_ora_vec = ome_ora$id
        res.ome = o.ea(ome_ora_vec, pathways,
                       uniq.count = unique.gene,
                       uniq.len = sapply(pathways, function(p) {p$gene.counts}))
        if (dim(res.meta$df)[1]==0 && dim(res.ome$df)[1]==0) { cat(glue("\nSkipping '{value_of_interest}' annotation-value; no significant enrichments.\n\n")); next }
        
        if (opt$pval_comb=="query") {
          if(print_internal_placemarks) cat("\n\n####################\nOverrepresenation Analysis on Both Omes\n\n")
          res.mo = o.ea(c(meta_ora_vec, ome_ora_vec), pathways,
                        uniq.count = unique.cmpd+unique.gene,
                        uniq.len = sapply(pathways, function(p) {p$cmpd.counts + p$gene.counts}))
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
                                  write_to_file = T, glue("{opt$output_prefix}_metabolome"))
      res.meta = q.ea(meta_qea$mat, meta_qea$cls, pathways, uniq.len = sapply(pathways, function(p) {p$cmpd.counts}))
      if (!multiomic && dim(res.meta$df)[1]==0) { cat(glue("\nSkipping '{value_of_interest}' annotation-value; no significant metabolite enrichments.\n\n")); next }
      
      if (multiomic) {
        if(print_internal_placemarks) cat(glue("\n\n####################\nQuantitative Enrichment Analysis on {opt$ome_type}\n\n"))
        ome_qea = gct.to.qea.input(ome_input, annot_of_interest, value_of_interest, annots = annots,
                                   write_to_file = T, prefix=glue("{opt$output_prefix}_{opt$ome_type}"))
        res.ome = q.ea(ome_qea$mat, ome_qea$cls, pathways, uniq.len = sapply(pathways, function(p) {p$gene.counts}))
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
          res.mo = q.ea(mo_qea$mat, mo_qea$cls, pathways, uniq.len = sapply(pathways, function(p) {p$cmpd.counts + p$gene.counts}))
        }
      }
      
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
          w.m <- unique.cmpd/(unique.cmpd+unique.gene)
          w.g <- unique.gene/(unique.cmpd+unique.gene)
        } else if (opt$pval_comb=="pvalp") {
          pw.c = sapply(all.paths, function(pathway) {pathways[[pathway]]$cmpd.counts})
          pw.g = sapply(all.paths, function(pathway) {pathways[[pathway]]$gene.counts})
          w.m = pw.c/(pw.c+pw.g)
          w.g = gene.w = pw.g/(pw.c+pw.g)
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
                            N.Entries = sapply(pathways, function(p) {p$cmpd.counts + p$gene.counts})[all.paths],
                            N.Hits = NA) %>% # initialize N.Hits column for order
          # add metabolite / gene specific hit breakdowns
          dplyr::mutate(N.Hits.Compounds = res.meta.df[ID,'N.Hits'] %>% ifelse(is.na(.), 0, .), # count 
                        N.Hits.Genes = res.ome.df[ID,'N.Hits'] %>% ifelse(is.na(.), 0, .),
                        N.Hits = N.Hits.Compounds+N.Hits.Genes,
                        Test.Type = opt$anal_type) %>%
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
                              res.meta$hits[all.paths], res.ome$hits[all.paths])
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
      ggplot(res.fin.df, aes(x = Impact.DC, y = neg.Log.P.Value, color = neg.Log.P.Value)) + 
        geom_point() +
        ylab("-log(P-Value)")+
        xlab("Impact (Degree Centrality)")
      # save to file
      fn = glue("{opt$output_prefix}_{make.names(annot_of_interest)}_{value_of_interest.name}_negLogFC_vs_Impact.png")
      ggsave(file.path(annot_dir, fn))
    }
    
    #### Enrichment Ratio ####
    enrichment_metric = "Impact.CC"
    plot_df = head(res.fin.df,opt$top_n_networks) %>%
      dplyr::arrange(neg.Log.P.Value)
    # make 'volcano' plot
    ggplot(plot_df, aes(x = !!sym(enrichment_metric),
                        y = factor(Pathway.Name, level = Pathway.Name), # sort by order of appearance
                        size = neg.Log.P.Value)) + 
      geom_point() + 
      ylab("Pathway")+
      xlab(enrichment_metric)
    # save to file
    fn = glue("{opt$output_prefix}_{make.names(annot_of_interest)}_{value_of_interest.name}_EnrichemntRatio.png")
    ggsave(file.path(annot_dir, fn))
    
    
    #### Network Graphs ####
    if (!is.null(graph_list)) {
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
    top_n_pathways = floor(max_hm_pathways/length(res.df.list)) # pathways to plot per subvalue
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
  hm <- Heatmap(heatmap_df, # plot normalized heatmap
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(is.na(sign_df[i, j])) {
                    grid.text("", x, y)
                  } else if(sign_df[i, j] <0.01) {
                    grid.text("***", x, y, gp = gpar(fontface = "bold", col = "white"))
                  } else if(sign_df[i, j] <0.02) {
                    grid.text("**", x, y, gp = gpar(fontface = "bold", col = "white"))
                  } else if(sign_df[i, j] <0.05) {
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





