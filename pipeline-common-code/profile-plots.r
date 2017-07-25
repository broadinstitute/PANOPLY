

library (reshape)
library (ggplot2)
library (RColorBrewer)


source ('preamble.r')

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


read.gct.data <- function (f) {
  d <- read.gct (f)
  return (d[,3:ncol(d)])
}



cls <- read.cls ( file.path (data.dir, bimodal.cls) )
data <- read.gct.data ( file.path (pre.dir, paste (type, '-ratio.gct', sep='')) )



## plot individual profiles
plot.profiles <- function (d.in, out.file, title) {
  d <- melt (d.in)
  colnames (d) <- c ('sample', 'ratio')
  s <- data.frame (sample=colnames (d.in), type=cls)
  d <- merge (d,s)
  g <- ggplot (aes (x=ratio, color=type), data=d) + geom_density(size=1) + xlab ('log ratio') + xlim(-8,4) +
    ggtitle (title) + facet_wrap (~sample, ncol=12) + theme_minimal() +
    # theme(panel.background=element_rect(fill='gray95')) + 
    scale_color_manual (name='Type', values=paired.colors (1, nlevels (d[,'type'])), 
                        labels=levels (d[,'type']))
  ggsave (out.file, width=28, height=21)
}

plot.profiles (data, paste (type, '-density-profiles.pdf', sep=''), 
               paste (toupper (type), ': Distribution of log ratios', sep=''))



## median ratio (SpectrumMill) profiles with mean ratio profiles

plot.comparative.profiles <- function (d.in, profile, plot=TRUE) {
  melt.cols <- c ('sample', 'ratio')
  d <- melt (d.in)
  colnames (d) <- melt.cols
  s <- data.frame (sample=colnames (d.in), type=cls)
  d.median <- merge (d,s)
  
  int.itraq <- read.gct.data ( file.path (pre.dir, paste (profile, '-itraq-intensity.gct', sep='')) )
  int.ref <- read.gct.data ( file.path (pre.dir, paste (profile, '-reference-intensity.gct', sep='')) )
  ratio.calc <- log2 (int.itraq / int.ref)
  
  d.r <- melt (ratio.calc)
  colnames (d.r) <- melt.cols
  d.mean <- merge (d.r, s)
  
  d.plot <- rbind (data.frame (d.median, method='median'), data.frame (d.mean, method='mean'))
  
  if (plot) {
    g <- ggplot (aes (x=ratio, color=interaction (type,method)), data=d.plot) + geom_density(size=1) + 
      xlab ('log ratio') + xlim(-8,4) + ggtitle (profile) + facet_wrap (~sample, ncol=12) + 
      # theme(panel.background=element_rect(fill='gray95')) +
      theme_minimal() +
      scale_color_manual (name='Type [Method]', values=paired.colors (2, nlevels (d.plot[,'type'])), 
                          labels = c (paste (levels (d.plot[,'type']), '[median]'),
                                      paste (levels (d.plot[,'type']), '[summed]')))
    
    ggsave (paste (profile, '-comparative-density-profiles.pdf', sep=''), width=28, height=21)
  }

  invisible (d.plot)
}

plot.comparative.profiles (data, type)


#
### For comparison with PNNL processing results -- uncomment when needed
### (only proteome data is available in crossTab.RData)
# ## median + mean + results from PNNL
# load ('../../all-tumor-analysis/bimodal-analysis-with-mi/pnnl/crossTab.RData')
# # set proper column names to match with other datasets
# cols.data <- colnames (data)
# cols <- sapply (colnames (crossTab), 
#                 function (x) { 
#                   # ensure replicate samples are properly handled
#                   s1 <- strsplit (x, split='\\.')[[1]]
#                   i <- ifelse (length(s1)>1, as.integer(s1[2]), 1)
#                   s <- grep (strsplit (s1[1], split='-')[[1]][3], cols.data)[i]
#                   return (s)
#                 })
# colnames (crossTab) <- cols.data [cols]
# pnnl <- melt (crossTab)[,2:3]
# colnames (pnnl) <- c ('sample', 'ratio')
# s <- data.frame (sample=colnames (data), type=cls)
# d.pnnl <- merge (pnnl, s)
# d.plot <- rbind (plot.comparative.profiles (data, 'proteome', plot=FALSE), data.frame (d.pnnl, method='pnnl'))
# # plot
# g <- ggplot (aes (x=ratio, color=interaction (type,method)), data=d.plot) + geom_density(size=1) + 
#   xlab ('log ratio') + xlim(-8,4) + ggtitle ('proteome') + facet_wrap (~sample, ncol=12) + 
#   theme(panel.background=element_rect(fill='gray95')) + 
#   scale_color_manual (name='Type [Method]', values=paired.colors (3), 
#                       labels=c('Bimodal [median]', 'Ignore [median]', 'Unimodal [median]',
#                                'Bimodal [summed]', 'Ignore [summed]', 'Unimodal [summed]',
#                                'Bimodal [pnnl]', 'Ignore [pnnl]', 'Unimodal [pnnl]'))
# 
# ggsave ('proteome-pnnl-comparative-density-profiles.pdf', width=28, height=21)
# 

