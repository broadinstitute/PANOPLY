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
                  ptmsigdb=T,                      ## if TRUE and PTMsigDB was used, separate heatmaps for the different 
                                                   ## PTMsigDB catagories witll be created (_KINASE, _PERT, _PATH, etc)
                  ser.meth='ARSA',                 ## Seriation method used to arrange the matrix. Only used if
                  cw=10,                           ## heatmap cellwidth 
                  ch=10,                           ## heatmap cellheight
                  
                  remove_prefix = FALSE,           ## Toggle to remove common prefixes for names (e.g. "HALLMARK_MITOTIC_SPINDLE" -> "MITOTIC_SPINDLE")
                  normalize_names = FALSE,        ## Toggle for normalizing pathway names (e.g. "MITOTIC_SPINDLE" -> "Mitotic Spindle")
                  
                  colors = c("#347db6", "#6babd0", "#afd3e6", '#FFFFFF', "#f6bda4", "#e48169", "#c33d3e") , ## default color-scale for heatmap (taken from Vega RdBu)
                  na.col = 'grey50',
                  # number_color = NULL,
                  
                  geneset_groups_file = NULL,      ## CSV with a mapping between genesets and groups
                  geneset_groups_colors = NULL,    ## named vector with colors corresponding to geneset groups
                  ...                              ## further argument passed onto pheatmap              
){
  library(pacman)
  p_load(RColorBrewer)
  p_load(pheatmap)
  p_load(seriation)
  p_load(glue)
  p_load(tidyverse)
  
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
    
    #### Identify Top Genesets below FDR threshold ####
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
    rdesc.filt <- rdesc[keep.idx,]
    
    #### Process Annotations ####
    
    ## add 'C' if column names are all numeric
    if(sum( is.na( suppressWarnings(as.numeric(colnames(mat.filt))) ) ) == 0) # suppress warning of as.numeric(), since we EXPECT this to be characters most of the time
      colnames(mat.filt) <- paste0('C', colnames(mat.filt))
    
    ## add significance star annotation
    anno.row <- matrix('', nrow=nrow(mat.filt), ncol=ncol(mat.filt), dimnames = dimnames(mat.filt))
    for(i in 1:ncol(mat.filt)) { # for each ssGSEA test
      idx.signif = which(fdr.filt[, i] < fdr.max) # determine which of the displayed features are significant
      anno.row[idx.signif, i] <- '*' # add a significance star to every significant test
    }
    
    ## get min / max values
    max.val = ceiling( max( abs(mat.filt), na.rm=T) )
    min.val = -max.val
    
    color.breaks = seq( min.val, max.val, length.out=100 )
    
    ## heatmap color
    # colors = c("#2166AC", "#5DA3CB", "#BBDAEA", "#F7F7F7", "#FAC9B0", "#E0775E", "#B2182B") # paul tol RdBu (khroma::color('BuRd')(7))
    # colors = c("#2166AC", "#3783BB", "#5DA3CB", "#92C5DE", "#F7F7F7", "#F4A582", "#E0775E", "#CA4841", "#B2182B") # paul tol RdBu middle colors removed
    # colors = c("#347db6", "#6babd0", "#afd3e6", "#e0e9ef", '#FFFFFF', "#f7e4d9", "#f6bda4", "#e48169", "#c33d3e") # taken from the vega red-blue color-palette
    # colors <- rev(RColorBrewer::brewer.pal(7, "RdBu")) # taken from RColorBrewer RdBu (reversed)
    # colors <- c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B") # taken from RColorBrewer RdBu (reversed), middle colors removed
    color.hm <- colorRampPalette(colors)(99) # taken from the vega red-blue color-palette
    
    ## reorder rows
    if('geneset_groups' %in% colnames(rdesc.filt)){ # if we added groups to our genesets
      ord.idx <- order(rdesc.filt$geneset_groups) # order by those cateogries
    } else { # otherwise attempt to cluster rows
      ord.idx = tryCatch({
        dist.row <- dist(mat.filt, method = 'euclidean') %>% # find distances
          seriate(dist.row, method = ser.meth) # compute order
        get_order(dist.row) # sort by that order
      }, error = function(e) { # if that fails
        message ("Could not cluster data-base signatures; heatmap rows will be left unsorted.\n")
        return (rownames(mat.filt))
      })
    }
    mat.filt <- mat.filt[ord.idx, , drop=F]
    anno.row <- anno.row[ord.idx, , drop=F]
    
    #### import process-categories, if available ####
    if('geneset_groups' %in% colnames(rdesc.filt)){
      rdesc.filt <- rdesc.filt[ord.idx, ]
      annotation_row <- matrix(rdesc.filt$geneset_groups, ncol=1, dimnames = list(rownames(rdesc.filt), c('Category')))
      annotation_row <- data.frame(annotation_row)
      annotation_colors = list(Category = as.vector(unique(geneset_groups_df[c("geneset_groups", "geneset_groups_colors")]))$geneset_groups_colors) # convert dataframe into named list of named vectors
      # annotation_colors <- list(Category=c('signaling'='skyblue2', 'immune'='coral2', 'development'='peachpuff', 
      #                                 'proliferation'='palegreen3', 'cellular component'='snow4', 'metabolic'='khaki', 
      #                                 'DNA damage'='darkmagenta', 'pathway'='tan3'))
      gaps_row <- cumsum(table(annotation_row$Category))
      if (is.hallmark) { # if is hallmark (legacy behavior)
        rownames(mat.filt) <- sub('^HALLMARK_', '',  rownames(mat.filt))
        rownames(annotation_row) <- sub('^HALLMARK_', '',  rownames(annotation_row))
      }
    } else {
      annotation_row <- NULL
      annotation_colors <- NULL
      gaps_row=NULL
    }   
    
    if (remove_prefix) {
      # remove common prefixes (e.g. "KEGG_MEDICUS_<pathway>" -> "<pathway>")
      prefixes = unique(gsub("^(.+?)_.+$", "\\1", rownames(mat.filt)))
      while ( dim(mat.filt)[1]>1 && length(prefixes)==1) { # if there's only one common prefix (ASSUMING WE HAVE MORE THAN ONE PATHWAY)
        rownames(mat.filt) <- sub(glue('^{prefixes}_'), '',  rownames(mat.filt)) # prune it from matrix
        if (!is.null(annotation_row)) rownames(annotation_row) <- sub(glue('^{prefixes}_'), '',  rownames(annotation_row)) # prune it annotations
        prefixes = unique(gsub("^(.+?)_.+$", "\\1", rownames(mat.filt))) # check for another common prefix
      }
    }
    if (normalize_names) {
      # replace underscores with spaces
      rownames(mat.filt) = gsub("_", " ", rownames(mat.filt))
      if (!is.null(annotation_row)) rownames(annotation_row) = gsub("_", " ", rownames(annotation_row))
      # replace underscores with spaces
      rownames(mat.filt) = str_to_title(rownames(mat.filt))
      if (!is.null(annotation_row)) rownames(annotation_row) = str_to_title(rownames(annotation_row))
    }
    
    # #### Creata Heatmap ####
    # try(pheatmap(mat.filt, 
    #              cluster_cols = F, 
    #              cluster_rows=F, 
    #              col=color.hm, 
    #              breaks = color.breaks, filename = fn.out, 
    #              display_numbers = anno.row, #number_color = number_color,
    #              na_col = na.col, cellwidth = cw, cellheight = ch, 
    #              annotation_row = annotation_row ,
    #              annotation_colors = annotation_colors,
    #              gaps_row = gaps_row,
    #              ...))
    
    
    
    #### Creata Heatmap ####
    try(pheatmap(mat.filt, 
                 cluster_cols = F, 
                 cluster_rows=F, 
                 col=color.hm, 
                 breaks = color.breaks, filename = fn.out, 
                 display_numbers = anno.row, #number_color = number_color,
                 na_col = na.col, cellwidth = cw, cellheight = ch, 
                 annotation_row = annotation_row ,
                 annotation_colors = annotation_colors,
                 gaps_row = gaps_row,
                 ...))
    
    
  }
  
  ################################################################
  ## PTMsigDB: separate heatmaps for different categories
  if(is.null(n.max)) n.max <- 'all'
  if(ptmsigdb){ ## split categories
    rid.type <- sub('^(.*?)-.*', '\\1', rid) %>% unique  
    for(rt in rid.type){
      fn.out=glue("heatmap_{rt}_max.fdr_{fdr.max}_n.max_{n.max}.pdf")  
      idx <- grep(glue("^{rt}"), rid)
      plothm(rdesc[idx, , drop=F], mat[idx, , drop=F], fdr.max, n.max, fn.out, cw, ch)
      fn.out=glue("heatmap_{rt}_max.fdr_{fdr.max}_n.max_{n.max}.png")  
      plothm(rdesc[idx, , drop=F], mat[idx, , drop=F], fdr.max, n.max, fn.out, cw, ch)
      
    }
  } else {
    fn.out=glue("heatmap_max.fdr_{fdr.max}_n.max_{n.max}.pdf")  
    tryCatch(plothm(rdesc, mat, fdr.max, n.max, fn.out, cw, ch),
             error = function(cond) {
               message("Unable to plot heatmaps, with the following Error:")
               message(paste(cond, "\n"))
             })
    fn.out=glue("heatmap_max.fdr_{fdr.max}_n.max_{n.max}.png")  
    tryCatch(plothm(rdesc, mat, fdr.max, n.max, fn.out, cw, ch),
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
    df = read.csv(geneset_groups_file)
    if ( dim(df)[2]!=2  ) {
      cat("WARNING: Geneset groups file is formatted incorrectly; must have two columns, with genesets in the first column and geneset-groupings in the second.") # warn user about geneset
      return(NULL)
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
  
  return( geneset_groups_df )
  
}

