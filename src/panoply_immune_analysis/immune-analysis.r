#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

source ('config.r')


if (!require("pacman")) install.packages ("pacman")
pacman::p_load (RColorBrewer, ComplexHeatmap, circlize, cluster, ape, tidyr, ggplot2, dplyr, devtools)
if (! require (estimate)) {
  # only available at R forge
  install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)
  library (estimate)
}
if (! require (xCell)) {
  # only available on GitHub
  devtools::install_github('dviraran/xCell')
  library (xCell)
}
if (! require (ImmuneSubtypeClassifier)) {
  # only available on GitHub
  devtools::install_github('Gibbsdavidl/ImmuneSubtypeClassifier')
  library (ImmuneSubtypeClassifier)
}

# 
# source ('/Volumes/prot_proteomics/LabMembers/manidr/R-utilities/io.r')
# rna.data.file <- 'rna.gct'
# immune.enrichment.subgroups <- 'groups-aggregate.csv'
# 

run.immune.analysis <- function (rna.data=rna.data.file, groups.file=immune.enrichment.subgroups,
                                 FDR=immune.enrichment.fdr,
                                 heatmap.width=10, heatmap.height=15) {
  
  ## read data
  rna <- parse.gctx (rna.data)
  
  
  ## estimate scores
  runEstimate <- function (ds, out.prefix) {
    # run estimate and save/return scores
    tmp.file <- sprintf ("tmp-%s.gct", Sys.getpid())
    out.file <- paste0 (out.prefix, '.gct')
    write.gct2 (data.frame (Name=rownames (ds), Description=rownames(ds), ds), tmp.file)
    estimateScore (tmp.file, out.file)
    
    scores.table <- parse.gctx (out.file)
    scores <- t (scores.table@mat)
    
    file.remove (tmp.file)
    file.remove (out.file)
    write.csv (scores, paste0 (out.prefix, '.csv'))
    return (scores)
  }
  rna.ES <- runEstimate (rna@mat, 'estimate-scores')
  
  
  ## xCell scores
  runXCell <- function (ds, out.prefix) {
    # run xCell and save/return scores
    scores.table <- xCellAnalysis (ds)
    scores <- t (scores.table)
    
    scores_outfile <- t (scores.table) %>% data.frame %>% tibble::rownames_to_column("Sample.ID")
    
    write.csv (scores_outfile, paste0 (out.prefix, '.csv'), row.names=FALSE)
    return (scores)
  }
  rna.XC <- runXCell (rna@mat, 'xcell-scores')
  
  
  ## PanTCGA-based Immune Subtype
  runImmuneSubtyping <- function (ds, out.prefix) {
    # run immune subtyping -- save/return subtype
    subtype <- callEnsemble (ds, geneids='symbol')
    colnames (subtype) <- c ('Sample.ID', 'Immune.Subtype', 
                             'Subtype1.WoundHealing',
                             'Subtype2.IFNGammaDominant',
                             'Subtype3.Inflammatory',
                             'Subtype4.LymphocyteDepleted',
                             'Subtype5.ImmunologicallyQuiet',
                             'Subtype6.TGFBetaDominant')
    write.csv (subtype, paste0 (out.prefix, '.csv'), row.names=FALSE)
    return (subtype)
  }
  rna.subtype <- runImmuneSubtyping (rna@mat, "immune-subtype")
  
  
  
  ## Enrichment Analysis
  # if subgroups file for enrichment is specified, run Fisher test
  if (! is.null (groups.file)) {
    # groups file specified for testing enrichments
    # (format similar to expt-design-file with Sample.ID and additional columns)
    # enrichement test will be run for each additional column
    
    subgroup.table <- read.csv (groups.file)
    rownames (subgroup.table) <- subgroup.table[,'Sample.ID']
    cls.list <- setdiff (colnames (subgroup.table), 'Sample.ID')
    
    rownames (rna.subtype) <- rna.subtype[,'Sample.ID']     
    subtype <- rna.subtype [subgroup.table [,'Sample.ID'], 'Immune.Subtype']    # match sample order 
    
    enrich <- NULL
    if (length (cls.list) > 0) {
      subtype.cls <- factor (subtype)
      for (cl in levels (subtype.cls)) {
        cl.cls <- ifelse (subtype.cls == cl, cl, "_rest_")
        cl.cls <- factor (cl.cls, levels=c(cl, '_rest_'))   # class of interest must be first level
        for (g in cls.list) {
          gl <- factor (subgroup.table [,g])
          for (gx in levels (gl)) {
            gx.cls <- ifelse (gl == gx, gx, "_other_")
            gx.cls <- factor (gx.cls, levels=c(gx, '_other_'))   # class of interest must be first level
            pval <- fisher.test (cl.cls, gx.cls, alternative='greater')$p.value
            enrich <- rbind (enrich, c (cl, g, gx, pval))
          }
        }
      }
    }
    
    # write out results
    colnames (enrich) <- c ('Immune.Subtype', 'group', 'subgroup', 'fisher.test.pvalue')
    enrich <- cbind (enrich, adj.pvalue=p.adjust (enrich[,'fisher.test.pvalue'], method='BH'))
    write.csv (enrich, "immune-subtype-enrichment.csv", row.names=FALSE)
    
    # extract statistically significant enrichments
    e <- enrich
    e.sig <- e [ e[,'fisher.test.pvalue'] < FDR, , drop=FALSE] 
    write.csv (e.sig, sprintf ('immune-subtype-enrichment-pval%.2f.csv', FDR), row.names=FALSE)
  }
  
  
  ###
  ### plots -- heatmap of xcell data, with annotation bars
  ###          xcell vs Estimate immune and stromal scores
  ###
  
  # heatmap
  results <- t (rna.XC [, -(1:3 + ncol (rna.XC))])
  if (!is.null (groups.file)) {
    samples <- as.character (subgroup.table[,'Sample.ID'])
    samples <- samples [ samples %in% colnames (results) ]
    annot <- subgroup.table [samples,]
  } else {
    samples <- colnames (results)
    annot <- data.frame (Sample.ID=samples)
  }
  rownames (annot) <- samples
  
  annot <- annot %>%
    mutate_at (names (.[, sapply (., function (x) length(unique(x)) <= 5)]), funs (factor(.))) %>%
    mutate (xCell.ImmuneScore = rna.XC[samples, 'ImmuneScore']) %>%
    mutate (xCell.StromaScore = rna.XC[samples, 'StromaScore']) %>% 
    mutate (xCell.MicroenvironmentScore = rna.XC[samples, 'MicroenvironmentScore']) %>%
    mutate (ESTIMATE.TumorPurity = rna.ES[samples, 'TumorPurity']) %>%
    mutate (Immune.Subtype = factor (rna.subtype[samples, 'Immune.Subtype'])) %>%
    select (-Sample.ID)
  annotation <- HeatmapAnnotation (df=annot, annotation_height = 0.5, annotation_width = 0.5,
                                   show_annotation_name=TRUE)
  heatmap <- Heatmap (results, top_annotation=annotation, 
                      row_names_gp = gpar (fontsize=8), column_names_gp = gpar(fontsize=8),
                      clustering_distance_columns='pearson', clustering_distance_rows='pearson',
                      width=heatmap.width, height=heatmap.height)
  pdf ('xcell-scores-heatmap.pdf', pointsize=8, 
       height=heatmap.height+10, width=heatmap.width+10)
  draw (heatmap)
  dev.off()
  
  # scatter plots
  plot.data <- data.frame (id=samples, xCell.ImmuneScore=rna.XC[samples,'ImmuneScore'], 
                           xCell.StromalScore=rna.XC[samples,'StromaScore'],
                           ESTIMATE.ImmuneScore=rna.ES[samples,'ImmuneScore'], 
                           ESTIMATE.StromalScore=rna.ES[samples,'StromalScore'])
  scatter.data <- gather (plot.data, key, value, -id) %>%
    separate (key, into=c('type', 'score'), sep='\\.') %>%
    spread (type, value)
  p <- ggplot (aes (x=ESTIMATE, y=xCell, group=score, color=score), data=scatter.data) + geom_point() +
    facet_wrap(~score, scales='free')
  ggsave ('xCell-vs-ESTIMATE-plots.pdf', width=8, height=4)
}


if (!interactive()) {
  # read parameters from config.r and run analysis
  rna.data <- rna.data.file
  groups.file <- ifelse (exists ("immune.enrichment.subgroups"), immune.enrichment.subgroups, NULL)
  FDR <- ifelse (exists ("immune.enrichment.fdr"), immune.enrichment.fdr, 0.05)  # default 0.05
  heatmap.width <- ifelse (exists ("immune.heatmap.width"), immune.heatmap.width, 10)
  heatmap.height <- ifelse (exists ("immune.heatmap.height"), immune.heatmap.height, 15)
  
  run.immune.analysis (rna.data=rna.data, groups.file=groups.file,
                       FDR=FDR, heatmap.width=heatmap.width, heatmap.height=heatmap.height)
}

