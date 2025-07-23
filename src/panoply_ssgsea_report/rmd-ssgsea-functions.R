#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
## ###########################################################################
##               PTM-SEA / ssGSEA volcano plots
gg_volc <- function(output.prefix, 
                    fdr.max = 0.05, ## max. FDR
                    n.max = 5, ## maximal number of signatures to label each side of the volcano plots
                    ncol=3,    ## number of plots per row
                    alpha=100,
                    #show=c('KINASE-PSP_mTOR/MTOR', 'KINASE-PSP_PKCA/PRKCA', 'KINASE-PSP_PKACA/PRKACA', 'KINASE-PSP_CDK2', 'KINASE-PSP_CDK1'),
                    show=NULL, ## signatures to show in all plots
                    ...  ## passed to ggsave
                    ){
    #fdr.max <- 0.05 
    #n.max <- 5 
    if(!require(pacman)) install.packages("pacman")
    p_load(ggpubr)
    p_load(ggrepel)
    p_load(dplyr)
  
    ## import ptm-sea results
    fn3 <- glue('{output.prefix}-combined.gct')
    ptm <- parse.gctx(fn3)
    mat <- ptm@mat
    cid <- ptm@cid
    rdesc <- ptm@rdesc
    rid <- ptm@rid
    
    
    ## loop over sample columns
    plot.list <- vector('list', length(cid))
    names(plot.list) <- cid
    
    for(s in cid){
      
      pval.tmp <- rdesc[ , glue("pvalue.{sub('^X','',make.names(s))}")] %>% as.numeric
      fdr.tmp <- rdesc[ , glue("fdr.pvalue.{sub('^X','',make.names(s))}")] %>% as.numeric
      ol.tmp <- rdesc[, glue("Signature.set.overlap.percent.{sub('^X','',make.names(s))}")] %>% as.numeric
      score.tmp <- mat[ , glue('{s}')] %>% as.numeric
      
      
      Enrichment.Score <- score.tmp
      P.Value <- -10*log(pval.tmp, 10)
      
      
      Significant <- rep('n', nrow(rdesc))
      sig.idx <- which(fdr.tmp < fdr.max) 
      if(length(sig.idx) > 0)
        Significant[sig.idx] <- 'y'
      
      col <- Significant
      col[col == 'y'] <- 'red'
      col[col == 'n'] <- 'grey'
      
      ## used for ranking
      Rank <- P.Value * Enrichment.Score
      
      data.plot <- data.frame(Enrichment.Score,
                              P.Value,
                              Significant,
                              Signature=rid,
                              Overlap.Percent=ol.tmp,
                              Rank=Rank,
                              col)
      if(!is.null(show)){
        highlight <- rep('n', nrow(data.plot))
        highlight[ which(data.plot$Signature %in% show )] <- 'y'
        data.plot <- data.frame(data.plot, highlight)
      }
      
      xlim=max(abs(data.plot$Enrichment.Score))
      ymax=max(abs(data.plot$P.Value))
      
      if(is.null(show)){
          data.sig.dn <- data.plot %>% 
            filter(Significant == 'y' & Enrichment.Score < 0) %>% 
            arrange(Rank) 
          if(nrow(data.sig.dn) > n.max)
            data.sig.dn <- data.sig.dn[1:(min( n.max, nrow(data.sig.dn))), ]
            #data.sig.dn <- data.sig.dn %>% slice(1:(min( n.max, nrow(data.sig.dn))))
          
          data.sig.up <- data.plot %>% 
            filter(Significant == 'y' & Enrichment.Score > 0) %>% 
            arrange(desc(Rank)) 
          if(nrow(data.sig.up) > n.max)
            data.sig.up <- data.sig.up[1:(min( n.max, nrow(data.sig.up))), ]
      } else {
        
        data.sig.up <- data.plot %>% 
          filter(highlight == 'y' & Enrichment.Score > 0) %>% 
          arrange(desc(Rank)) 
        data.sig.dn <- data.plot %>% 
          filter(highlight == 'y' & Enrichment.Score < 0) %>% 
          arrange(Rank) 
      }
      #  data.sig.up <- data.sig.up %>% slice(1:(min( n.max, nrow(data.sig.up))))
      
      
      p <- ggplot(data.plot, aes(Enrichment.Score, P.Value ) ) + 
        geom_point(aes(size=Overlap.Percent, colour=Significant)) + 
        xlim(-xlim, xlim) +
        geom_vline(xintercept=0, linetype=3) +
        labs(title=glue('{s}'), y='-10log10(P-value)') +
        scale_colour_manual(name="", values = c("y"=my.col2rgb("darkred", alpha = alpha), 
                                                "n"=my.col2rgb("grey", alpha = alpha))) +
        geom_label_repel( aes(label=Signature),
                          #nudge_x=data.sig.dn$Enrichment.Score,
                          nudge_x=-5,
                          #nudge_y=seq(max(data.plot$P.Value), min(data.plot$P.Value), length.out = nrow(data.sig.dn) )[order(data.sig.dn$P.Value, decreasing = T)],
                          nudge_y=-5,
                          direction='y',
                          data=data.sig.dn,
                          force=1, size=2) +
        
        geom_label_repel( aes(label=Signature),
                          #nudge_x=data.sig.up$Enrichment.Score,
                          nudge_x=5,
                          nudge_y=-5,
                          direction='y',
                          data=data.sig.up,
                          force=1, size=2)
      
      p <- p + theme_bw() + theme(plot.title=element_text(hjust=0.5))
      
      
      plot.list[[s]] <- p
    }
    
    ## plot and save
    if(length(plot.list) > 9){
      
      pdf(glue("volcano_{output.prefix}.pdf"), ...)
      for(s in names(plot.list)){
        plot(plot.list[[s]])
      }
      dev.off()
      
    } else {
      pp <- ggarrange(plotlist= plot.list, common.legend=T, ncol = ncol)
      ggsave(glue("volcano_{output.prefix}.pdf"),  device='pdf', ...)
    }
    #save(pp, file='debug.RData')
    #annotate_figure(pp, top=glue("class vector: {output.prefix}"))
    
    
   # return(data.plot)
   return(0)
}

## ######################################################
## pathway heatmap
## - developed for GSEA/PTM-SEA on NMF results
## 
pw_hm <- function(output.prefix, 
                  fdr.max = 0.05,                  ## max. FDR
                  n.max = NULL,                    ## maximal number of signatures to label each side of the volcano plots
                                                   ## ignired if set to NULL
                  split.by.prefix = NULL,          ## if TRUE, separate heatmaps will be created for pathways with unique prefixes "<prefix>-<pathway_name>"
                                                   ## (e.g. 'KINASE-PSP_CDC7' will be grouped with other 'KINASE' pathways)
                                                   ## generalized replacement for the `ptmsigdb` parameter; will automatically be set to TRUE if ptmsigdb is detected
                  cluster.rows=TRUE,               ## Toggle for clustering rows using distance matrix, with ser.meth as seriation method
                  ser.meth='ARSA',                 ## Seriation method used to (attempt to) arrange the matrix.
                  cw=10,                           ## heatmap cellwidth 
                  ch=10,                           ## heatmap cellheight
                  
                  remove_prefix = FALSE,           ## Toggle to remove common prefixes for names (e.g. "HALLMARK_MITOTIC_SPINDLE" -> "MITOTIC_SPINDLE")
                  normalize_names = FALSE,         ## Toggle for normalizing pathway names (e.g. "MITOTIC_SPINDLE" -> "Mitotic Spindle")
                  
                  colors = c("#347db6", "#6babd0", "#afd3e6", '#FFFFFF', "#f6bda4", "#e48169", "#c33d3e") , ## default color-scale for heatmap (taken from Vega RdBu)
                                                   ## Use "legacy" for blue -> orange colors scheme
                  na.col = 'grey50',
                  # number_color = NULL,
                  
                  geneset_groups_file = NULL,      ## CSV with a mapping between genesets and groups
                  geneset_groups_colors = NULL,    ## named vector with colors corresponding to geneset groups
                  ...                              ## further argument passed onto pheatmap              
){
  library(pacman)
  p_load(RColorBrewer)
  # p_load(pheatmap)
  p_load(seriation)
  p_load(glue)
  p_load(tidyverse)
  # new dependencies
  p_load(ComplexHeatmap)
  p_load(circlize)
  p_load(grid)
  p_load(gridExtra)
  
  #### import ssGSEA results ####
  if(class(output.prefix) == 'character'){
    fn3 <- glue('{output.prefix}-combined.gct')
    gct <- parse.gctx(fn3)
  } else {
    gct <- output.prefix
  }
  mat <- gct@mat
  cid <- gct@cid
  rdesc <- gct@rdesc
  rid <- gct@rid
  
  #### import geneset groups ####
  is.hallmark = all(grepl("^HALLMARK_", rid))
  if (is.null(geneset_groups_file) && is.hallmark) { use_legacy_hallmarks = TRUE } else { use_legacy_hallmarks = FALSE }
  geneset_groups_df = process_geneset_groups(geneset_groups_file = geneset_groups_file,
                                             geneset_groups_colors = geneset_groups_colors,
                                             use_legacy_hallmarks = use_legacy_hallmarks)
  
  #### if we have geneset_groups, add that annotation to rdesc ####
  if( !is.null(geneset_groups_df) ){
    
    if ( length(intersect(rid, rownames(geneset_groups_df))) == 0 ) {
      cat("WARNING: No overlap found between provided genesets in ssGSEA results and provided geneset groups. Genesets will be plotted without groups.") # warn user about geneset
    } else {
      # rbind geneset annotations to rdesc
      rdesc = left_join(tibble::rownames_to_column(rdesc),
                        tibble::rownames_to_column(geneset_groups_df),
                        by = 'rowname') %>%
        tibble::column_to_rownames('rowname')
    }

  }
  
  #########################
  ## helper function
  plothm <- function(rdesc, mat, fdr.max, n.max, fn.out, cw, ch){
    
    #### Data Filtering ####
    ## Identify top genesets below FDR thresh
    fdr <- rdesc[,grep('^fdr.pvalue', colnames(rdesc)), drop=F] 
    keep.idx.list <- lapply(1:ncol(fdr), function(i, fdr, mat){
     # cat(i)
      f=fdr[, i] ## fdr
      s=mat[, i] ## score
      idx=which(f < fdr.max)
      if(length(idx) > 0 & !is.null(n.max))
        idx=idx[order(abs(s[idx]), decreasing = T)[1:min(n.max, length(idx))]]
      rownames(mat)[idx]
    }, fdr, mat )
    keep.idx <- unique(unlist(keep.idx.list))
    if (length(keep.idx)==0) { cat(glue("## WARNING: No features found with FDR below selected threshold ({fdr.max})")); if (exists("rt")) {cat(glue(" for ID-Type {rt}\n"))}; return() } # stop function if matrix is 0x0; color-gen will error
    
    ## filter dataframes to those genesets
    fdr.filt <- fdr[keep.idx, , drop=F ]
    mat.filt <- mat[keep.idx, , drop=F]
    rdesc.filt <- rdesc[keep.idx, , drop=F]
    
    ## add 'C' if column names are all numeric
    if(sum( is.na( suppressWarnings(as.numeric(colnames(mat.filt))) ) ) == 0) # suppress warning of as.numeric(), since we EXPECT this to be characters most of the time
      colnames(mat.filt) <- paste0('C', colnames(mat.filt))
    
    
    #### Heatmap Colors ####
    
    ## get min / max values
    max.val = ceiling( max( abs(mat.filt), na.rm=T) )
    min.val = -max.val
    
    
    ## color scale
    # colors = c("#2166AC", "#5DA3CB", "#BBDAEA", "#F7F7F7", "#FAC9B0", "#E0775E", "#B2182B") # paul tol RdBu (khroma::color('BuRd')(7))
    # colors = c("#2166AC", "#3783BB", "#5DA3CB", "#92C5DE", "#F7F7F7", "#F4A582", "#E0775E", "#CA4841", "#B2182B") # paul tol RdBu middle colors removed
    # colors = c("#347db6", "#6babd0", "#afd3e6", "#e0e9ef", '#FFFFFF', "#f7e4d9", "#f6bda4", "#e48169", "#c33d3e") # taken from the vega red-blue color-palette
    # colors <- rev(RColorBrewer::brewer.pal(7, "RdBu")) # taken from RColorBrewer RdBu (reversed)
    # colors <- c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B") # taken from RColorBrewer RdBu (reversed), middle colors removed
    # color.hm <- colorRampPalette(colors)(99) # taken from the vega red-blue color-palette
    
    if (length(colors)==1 && colors == "legacy") colors = c('cyan','darkblue', 'grey90', 'orange', 'yellow')
    
    ## color mapping
    color.breaks = seq( min.val, max.val, length.out=length(colors) )
    col_fun <- circlize::colorRamp2(color.breaks, colors)
    
    
    #### Row Reordering ####
    if(!cluster.rows){ # if we turn off row-clustering
      ord.idx <- 1:dim(mat.filt)[1] # retain original order
    } else { # otherwise attempt to cluster rows
      ord.idx = tryCatch({ # always attempt to order rows; row_split will take care of groupings
        dist.row <- dist(mat.filt, method = 'euclidean') %>% # find distances
          seriate(method = ser.meth) # compute order
        get_order(dist.row) # sort by that order
      }, error = function(e) { # if that fails
        message ("Could not cluster data-base signatures; heatmap rows will be left unsorted.\n")
        return (rownames(mat.filt))
      })
    }
    mat.filt <- mat.filt[ord.idx, , drop=F]
    fdr.filt <- fdr.filt[ord.idx, , drop=F]
    rdesc.filt <- rdesc.filt[ord.idx, , drop=F]
    
    
    #### Rowname Formatting ####
    if (remove_prefix) {
      # remove common prefixes (e.g. "KEGG_MEDICUS_<pathway>" -> "<pathway>")
      prefixes = unique(gsub("^(.+?)_.+$", "\\1", rownames(mat.filt)))
      while ( dim(mat.filt)[1]>1 && length(prefixes)==1) { # if there's only one common prefix (ASSUMING WE HAVE MORE THAN ONE PATHWAY)
        rownames(mat.filt) <- sub(glue('^{prefixes}_'), '',  rownames(mat.filt)) # prune it from matrix
        prefixes = unique(gsub("^(.+?)_.+$", "\\1", rownames(mat.filt))) # check for another common prefix
      }
    }
    if (is.hallmark) rownames(mat.filt) <- sub('^HALLMARK_', '',  rownames(mat.filt)) # if is hallmark, also prune prefix (legacy behavior)
    
    if (normalize_names) {
      # replace underscores with spaces
      rownames(mat.filt) = gsub("_", " ", rownames(mat.filt))
      # replace underscores with spaces
      rownames(mat.filt) = str_to_title(rownames(mat.filt))
    }
    
    
    
    #### Pathway Groupings (if provided) ####
    row_ha <- NULL
    row_split <- NULL
    if ('geneset_groups' %in% colnames(rdesc.filt)) {
      row_split <- rdesc.filt$geneset_groups
      # Optionally, add color annotation
      group_colors <- unique(geneset_groups_df$geneset_groups_colors)
      names(group_colors) <- unique(geneset_groups_df$geneset_groups)
      row_ha <- rowAnnotation(
        Category = row_split,
        col = list(Category = group_colors)
      )
    }
    
    ## add significance star annotation
    cell_fun <- function(j, i, x, y, width, height, fill) {
      if (fdr.filt[i, j] < fdr.max) {
        gb = textGrob("*")
        gb_w = convertWidth(grobWidth(gb), "mm")
        gb_h = convertHeight(grobHeight(gb), "mm")
        grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
      }
    }

    #### Create Heatmap #### 
    # Draw the heatmap
    ht <- Heatmap(
      mat.filt,
      name = "NES",
      col = col_fun,
      cluster_rows = FALSE, # using distance matrix and seriation method earlier in the code
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      na_col = na.col,
      row_split = row_split,
      row_title=NULL,
      left_annotation = row_ha,
      cell_fun = cell_fun,
      width = unit(ncol(mat.filt) * cw*1.5, "point"),
      height = unit(nrow(mat.filt) * ch*1.5, "point")
    )
    # Create heatmap (without legends)
    ht_fig <- draw(ht, 
                   show_heatmap_legend = FALSE, 
                   show_annotation_legend = FALSE)
    ht_grob = grid.grabExpr(draw(ht,
                                 show_heatmap_legend = FALSE, 
                                 show_annotation_legend = FALSE))
    # # Wrap in left-aligned viewport
    # ht_left_justified <- grobTree(
    #   ht_grob,
    #   vp = viewport(x = unit(0, "npc"), just = "left")
    # )
    
    # Create legend(s) figure
    base_legend = ComplexHeatmap::Legend(title = "NES", col_fun = col_fun) # get NES legend (used in both legends)
    min_legend_width = ComplexHeatmap::width.Legends(ComplexHeatmap::Legend(title = "NES", col_fun = col_fun))
    if (!is.null(row_ha)) {
      lgd_fig <- packLegend(base_legend,
                            ComplexHeatmap::Legend(title = "Category", at = names(group_colors), legend_gp = gpar(fill = group_colors)),
                            direction = "vertical")
      max_legend_text = max_text_width( c("Category", names(group_colors)) ) # get max text_width of groups
    } else {
      lgd_fig <- base_legend
      max_legend_text = max_text_width( "Category" ) # use the word "Category" as placeholder for max text_width
    }
    lgd_grob = grid.grabExpr(draw(lgd_fig, x=unit(0.95, "npc"), just="right")) # align legend as far to the right as possible
    
    #### Plot Heatmap & Legend #### 
    for (ext in c('.pdf', '.png')) { ## create pdf and png
      padding = 1 # padding to add to the width of each figure
      hm_width = convertX(unit(ncol(mat.filt) * cw*1.5, "point"), 'inches', valueOnly = TRUE) + # heatmap width
        convertX(ComplexHeatmap::max_text_width(rownames(mat.filt)), 'inches', valueOnly = TRUE) # rowname label width, +1 for padding
      lgd_width = convertX(min_legend_width, 'inches', valueOnly = TRUE) + # add minimum legend width
        convertX(max_legend_text, 'inches', valueOnly = TRUE) # add max legend text-length
      width = hm_width+lgd_width + padding # calculate total width + padding
        
      height = convertX(unit(nrow(mat.filt) * ch*1.5, "point"), 'inches', valueOnly = TRUE) + # heatmap height
        convertX(ComplexHeatmap::max_text_width(colnames(mat.filt)),'inches', valueOnly = TRUE) + padding # column-name text, + 1 for padding
      if (ext=='.pdf') pdf(paste0(fn.out,ext), width = width, height = height)
      if (ext=='.png') png(paste0(fn.out,ext), width = width, height = height, units = 'in', res=300)
      
      ## plot heatmap and legend side-by-side
      grid.arrange(ht_grob, lgd_grob, ncol=2, widths=c(hm_width,lgd_width))
      
      dev.off()
    }
    
  }
  
  ################################################################
  ## Plot Heatmaps
  
  ## check if dataset is PTM-SEA
  is_ptmsea = mean(grepl("^(PERT-PSP)|(PERT-P100-PRM)|(PERT-P100-DIA)|(PATH-WP)|(PATH-NP)|(KINASE-PSP)|(DISEASE-PSP)_.+?$", rid)) > 0.9  # Heuristic: at least 90% of rids fit the PTM-SEA prefixes. realistically, should be 100%, but adding flexibility in case the database changes slightly
  ## split into multiple heatmaps based on pathway prefixes, if we have ptmsea
  if (is.null(split.by.prefix) && is_ptmsea) split.by.prefix=TRUE # if we didn't set split.by.prefix explicitly, use heuristic to determine if we wanna split
  if (!is.null(split.by.prefix) && split.by.prefix) { # if split.by.prefix
    rid.type <- sub('^(.*?)-.*', '\\1', rid) %>% unique # split rid into groups according to prefixes
  } else { rid.type = "" }
  
  ## if na.max is NULL, use 'all' for n.max label in filename
  if(is.null(n.max)) { n.max.fn <- 'all' } else { n.max.fn = n.max }
  
  ## plot heatmaps
  for (rt in rid.type) {
    idx <- grep(glue("^{rt}"), rid) # filter dataset to pathway-grouping
    rt.fn = ifelse(rt=="", rt, paste0('_',rt)) # label for filename
    
    ## plot heatmap as PDF and PNG
    fn.out=glue("heatmap{rt.fn}_max.fdr_{fdr.max}_n.max_{n.max.fn}")
    tryCatch(plothm(rdesc[idx, , drop=F], mat[idx, , drop=F],
                    fdr.max, n.max,
                    fn.out,
                    cw, ch),
             error = function(cond) {
               message("Unable to plot heatmaps, with the following Error:")
               message(paste(cond, "\n"))
             })
  }
  
}




process_geneset_groups <- function(geneset_groups_file=NULL,
                                   geneset_groups_colors = NULL, 
                                   use_legacy_hallmarks = FALSE) {
  # default color-palette taken from Paul Tol light
  if (is.null(geneset_groups_colors)) geneset_groups_colors = c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#BBCC33", "#AAAA00", "#DDDDDD") # (optionally named) vector with color-palette for groups; 
  #### Create geneset_groups Object ####
  # if we have no groups_file AND we aren't using the legacy hallmarks colors
  if (is.null(geneset_groups_file) && !use_legacy_hallmarks) { geneset_groups=NULL } # set both geneset objects to NULL
  # if we have a groups file
  if (!is.null(geneset_groups_file)) {
    # use that groups.file to set the groups
    df = tryCatch({ read.csv(geneset_groups_file) }, error = function(e) { # if that fails
            warning (glue("Could not read in geneset-groupings file '{geneset_groups_file}' as a CSV. Genesets will be left ungrouped.\n"))
            return (NULL)
          })
    if (is.null(df)) return(NULL) # return null if geneset groupings file is missing
                  
    if ( dim(df)[2]!=2  ) {
      cat("WARNING: Geneset groups file is formatted incorrectly; must have two columns, with genesets in the first column and geneset-groupings in the second.") # warn user about geneset
      return(NULL) # return null if geneset groupings file is malformed
    }
    # if all checks were passed, convert to named vector
    geneset_groups = tibble::deframe(df)
  }
  
  # if we're using legacy hallmarks colors
  if (use_legacy_hallmarks) {
    ###################################
    geneset_groups <- c(
      HALLMARK_TNFA_SIGNALING_VIA_NFKB='signaling',
      HALLMARK_HYPOXIA='pathway',
      HALLMARK_CHOLESTEROL_HOMEOSTASIS='metabolic',
      HALLMARK_MITOTIC_SPINDLE='proliferation',
      HALLMARK_WNT_BETA_CATENIN_SIGNALING='signaling',
      HALLMARK_TGF_BETA_SIGNALING='signaling',
      HALLMARK_IL6_JAK_STAT3_SIGNALING='immune',
      HALLMARK_DNA_REPAIR='DNA damage',
      HALLMARK_G2M_CHECKPOINT='proliferation',
      HALLMARK_APOPTOSIS='pathway',
      HALLMARK_NOTCH_SIGNALING='signaling',
      HALLMARK_ADIPOGENESIS='development',
      HALLMARK_ESTROGEN_RESPONSE_EARLY='signaling',
      HALLMARK_ESTROGEN_RESPONSE_LATE='signaling',
      HALLMARK_ANDROGEN_RESPONSE='signaling',
      HALLMARK_MYOGENESIS='development',
      HALLMARK_PROTEIN_SECRETION='pathway',
      HALLMARK_INTERFERON_ALPHA_RESPONSE='immune',
      HALLMARK_INTERFERON_GAMMA_RESPONSE='immune',
      HALLMARK_APICAL_JUNCTION='cellular component',
      HALLMARK_APICAL_SURFACE='cellular component',
      HALLMARK_HEDGEHOG_SIGNALING='signaling',
      HALLMARK_COMPLEMENT='immune',
      HALLMARK_UNFOLDED_PROTEIN_RESPONSE='pathway',
      HALLMARK_PI3K_AKT_MTOR_SIGNALING='signaling',
      HALLMARK_MTORC1_SIGNALING='signaling',
      HALLMARK_E2F_TARGETS='proliferation',
      HALLMARK_MYC_TARGETS_V1='proliferation',
      HALLMARK_MYC_TARGETS_V2='proliferation',
      HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION='development',
      HALLMARK_INFLAMMATORY_RESPONSE='immune',
      HALLMARK_XENOBIOTIC_METABOLISM='metabolic',
      HALLMARK_FATTY_ACID_METABOLISM='metabolic',
      HALLMARK_OXIDATIVE_PHOSPHORYLATION='metabolic',
      HALLMARK_GLYCOLYSIS='metabolic',
      HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY='pathway',
      HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY='pathway',
      HALLMARK_P53_PATHWAY='proliferation',
      HALLMARK_UV_RESPONSE_UP='DNA damage',
      HALLMARK_UV_RESPONSE_DN='DNA damage',
      HALLMARK_ANGIOGENESIS='development',
      HALLMARK_HEME_METABOLISM='metabolic',
      HALLMARK_COAGULATION='immune',
      HALLMARK_IL2_STAT5_SIGNALING='signaling',
      HALLMARK_BILE_ACID_METABOLISM='metabolic',
      HALLMARK_PEROXISOME='cellular component',
      HALLMARK_ALLOGRAFT_REJECTION='immune',
      HALLMARK_SPERMATOGENESIS='development',
      HALLMARK_KRAS_SIGNALING_UP='signaling',
      HALLMARK_KRAS_SIGNALING_DN='signaling',
      HALLMARK_PANCREAS_BETA_CELLS='development'
    )
    ###################################
    
    geneset_groups_colors <- c('signaling'='skyblue2', 'immune'='coral2', 'development'='peachpuff', 
                               'proliferation'='palegreen3', 'cellular component'='snow4', 'metabolic'='khaki', 
                               'DNA damage'='darkmagenta', 'pathway'='tan3')
  }
  
  ## Return NULL if no genesets
  if (is.null(geneset_groups)) return(geneset_groups) # return NULL if no annotations exist
  
  
  ## Otherwise, Wrangle genesets into a Dataframe
  geneset_groups_df = data.frame(geneset_groups = geneset_groups) #otherwise create dataframe with geneset_groups
  
  if ( length(intersect(names(geneset_groups_colors), geneset_groups)) == 0 ) { # if we tried to set colors, but none of the color-names overlap with our set
    if ( !is.null(names(geneset_groups_colors)) )  cat("WARNING: Geneset groups and geneset color-scales share no overlapping values! Please check inputs; colors will be assigned arbitrarily.") # warn user about geneset
    if ( length(unique(geneset_groups)) > length(geneset_groups_colors) ) stop(paste0("Color-palette has too few values (",
                                                                                      length(geneset_groups_colors), ") ",
                                                                                      "for number of genesets (",length(unique(geneset_groups)),")"))
    # subset unnamed colors to length of geneset_groups & and name vector
    geneset_groups_colors = geneset_groups_colors[1:length(unique(geneset_groups))]
    names(geneset_groups_colors) = unique(geneset_groups)
  }
  geneset_groups_df = mutate(geneset_groups_df,
                             geneset_groups_colors = geneset_groups_colors[geneset_groups])
  rownames(geneset_groups_df) = names(geneset_groups) # reassign rownames
  
  return( geneset_groups_df )
  
}

