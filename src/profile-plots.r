#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

library (reshape)
library (ggplot2)
library (RColorBrewer)

source ('config.r')


paired.colors <- function (sets, types) {
  # types is 1, 2, or 3
  
  reds <- brewer.pal (7, 'Reds')
  greens <- brewer.pal (7, 'Greens')
  blues <- brewer.pal (7, 'Blues')
  
  colors <- NULL
  for (i in seq (7, (7-2*(sets-1)), -2)) {
    colors <- c (colors, reds[i])
    if (types==3) colors <- c (colors, greens[i], blues[i])
    else if (types==2) colors <- c (colors, greens[i])
  }
  
  return (colors)
}



## plot individual profiles
plot.profiles <- function (d.in, out.file, title, sampleinfo, width=28, height=21) {
  d <- melt (d.in)
  colnames (d) <- c ('sample', 'ratio')
  colnames (sampleinfo) <- c ('sample', 'type')
  
  s <- sampleinfo
  d <- merge (d,s)
  g <- ggplot (aes (x=ratio, color=type), data=d) + geom_density(size=1) + xlab ('log ratio') + # xlim(-8,4) +
    ggtitle (title) + facet_wrap (~sample, ncol=12) + theme_minimal() +
    # theme(panel.background=element_rect(fill='gray95')) + 
    scale_color_manual (name='Type', values=paired.colors (1, nlevels (d[,'type'])), 
                        labels=levels (d[,'type']))
  ggsave (out.file, width=width, height=height)
}



data <- parse.gctx ( file.path (pre.dir, paste (type, '-ratio.gct', sep='')) )
if (qc.col %in% colnames (data@cdesc)) { 
  sinfo <- data@cdesc [, c ('Sample.ID', qc.col)]
} else {
  if (!is.null (sampleQC.cls)) sinfo <- data.frame (data@cid, read.cls ( file.path (data.dir, sampleQC.cls) ))
  else sinfo <- data.frame (data@cid, rep (qc.pass.label))
}
plot.profiles (data@mat, paste (type, '-density-profiles.pdf', sep=''), 
               paste (toupper (type), ': Distribution of log ratios', sep=''), 
               sinfo)

