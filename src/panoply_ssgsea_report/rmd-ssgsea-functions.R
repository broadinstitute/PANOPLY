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
                  fdr.max = 0.05,  ## max. FDR
                  n.max = NULL,    ## maximal number of signatures to label each side of the volcano plots
                                   ## ignired if set to NULL
                  ptmsigdb=T,      ## if TRUE and PTMsigDB was used, separate heatmaps for the different 
                                   ## PTMsigDB catagories witll be created (_KINASE, _PERT, _PATH, etc)
                  ser.meth='ARSA', ## Seriation method used to arrange the matrix. Only used if
                  cw=10,           ## heatmap cellwidth 
                  ch=10,           ## heatmap cellheight
                  ...              ## further argument passed onto pheatmap              
){
  library(pacman)
  p_load(RColorBrewer)
  p_load(pheatmap)
  p_load(seriation)
  
  ## import ptm-sea results
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
  
  ###################################
  hallmark_process_category <- c(
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
  ## check if hallmark database was used  
  if(sum(rid %in% names(hallmark_process_category)) > 0){
    process_category <- rep('', length(rid))
    names(process_category) <- rid
    
    keep.idx <- which(names(hallmark_process_category) %in% names(process_category))
    hallmark_process_category <- hallmark_process_category[keep.idx]
    
    process_category[names(hallmark_process_category)] <- hallmark_process_category
    
    rdesc <- data.frame(rdesc, process_category)
  }
  
  #########################
  ## helper function
  plothm <- function(rdesc, mat, fdr.max, n.max, fn.out, cw, ch){
    
    ## fdr & score
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
    
    fdr.filt <- fdr[keep.idx, , drop=F ]
    mat.filt <- mat[keep.idx, , drop=F]
    rdesc.filt <- rdesc[keep.idx,]
    
    ## add 'C' if column names are all numeric
    if(sum( is.na(as.numeric(colnames(mat.filt))) ) == 0)
      colnames(mat.filt) <- paste0('C', colnames(mat.filt))
    
    anno.row <- matrix('', nrow=nrow(mat.filt), ncol=ncol(mat.filt), dimnames = dimnames(mat.filt))
    for(i in 1:ncol(mat.filt))
      anno.row[keep.idx.list[[i]], i] <- '*'
    
    dist.row <- dist(mat.filt, method = 'euclidean')

    max.val = ceiling( max( abs(mat.filt), na.rm=T) )
    min.val = -max.val
    
    color.breaks = seq( min.val, max.val, length.out=100 )
    
    ## heatmap color
    color.hm <- colorRampPalette(c('cyan','darkblue', 'grey90', 'orange', 'yellow'))(99)
  
    if('process_category' %in% colnames(rdesc.filt)){
      ord.idx <- order(rdesc.filt$process_category)
    } else {
      dist.row <- seriate(dist.row, method = ser.meth)
      ord.idx <- get_order(dist.row)
    }
    mat.filt <- mat.filt[ord.idx, , drop=F]
    anno.row <- anno.row[ord.idx, , drop=F]
    
    if('process_category' %in% colnames(rdesc.filt)){
      rdesc.filt <- rdesc.filt[ord.idx, ]
      annotation_row <- matrix(rdesc.filt$process_category, ncol=1, dimnames = list(rownames(rdesc.filt), c('cat')))
      annotation_row <- data.frame(annotation_row)
      annotation_colors <- list(cat=c('signaling'='skyblue2', 'immune'='coral2', 'development'='peachpuff', 
                                      'proliferation'='palegreen3', 'cellular component'='snow4', 'metabolic'='khaki', 
                                      'DNA damage'='darkmagenta', 'pathway'='tan3'))
      gaps_row <- cumsum(table(annotation_row$cat))
      rownames(mat.filt) <- sub('^HALLMARK_', '',  rownames(mat.filt))
      rownames(annotation_row) <- sub('^HALLMARK_', '',  rownames(annotation_row))
        
    } else {
      annotation_row <- NA
      annotation_colors <- NA
      gaps_row=NULL
    }   
    
    ## creata heatmap
    try(pheatmap(mat.filt, 
                 cluster_cols = F, 
                 cluster_rows=F, 
                 col=color.hm, 
                 breaks = color.breaks, filename = fn.out, 
                 display_numbers = anno.row, 
                 na_col = 'white', cellwidth = cw, cellheight = ch, 
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
      plothm(rdesc[idx, ], mat[idx, ], fdr.max, n.max, fn.out, cw, ch)
      fn.out=glue("heatmap_{rt}_max.fdr_{fdr.max}_n.max_{n.max}.png")  
      plothm(rdesc[idx, ], mat[idx, ], fdr.max, n.max, fn.out, cw, ch)
      
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
