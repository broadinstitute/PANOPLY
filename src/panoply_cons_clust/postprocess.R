#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

source ('config.r')

best.clus <- read.table ( sprintf ("%s.bestclus.txt", type), header=TRUE)
best.clus <- best.clus [,-3]
colnames (best.clus) <- c ('Sample.ID', 'Cluster')

write.csv (best.clus, sprintf ("%s-bestclus.csv", type), row.names=FALSE, quote=FALSE)

# add groups file to config.r for running association analysis
write ( sprintf ("assoc.subgroups <- '%s-bestclus.csv'", type), file="config.r", append=TRUE )


## Enrichment Analysis
# if subgroups file for enrichment is specified, run Fisher test
if (exists ("cluster.enrichment.subgroups")) {
  # groups file specified for testing enrichments
  # (format similar to expt-design-file with Sample.ID and additional columns)
  # enrichement test will be run for each additional column
  
  subgroup.table <- read.csv (cluster.enrichment.subgroups)
  rownames (subgroup.table) <- subgroup.table [, 'Sample.ID']
  cls.list <- setdiff (colnames (subgroup.table), 'Sample.ID')
  
  gct.file <- file.path (norm.dir, paste (master.prefix, '.gct', sep=''))
  ds <- parse.gctx (gct.file)
  sample.order <- ds@cid
  subgroup.table <- subgroup.table [ sample.order, ]     # match sample order with cluster cls vector
  
  enrich <- NULL
  if (length (cls.list) > 0) {
    cluster.cls <- factor (best.clus [,'Cluster'])
    for (cl in levels (cluster.cls)) {
      cl.cls <- ifelse (cluster.cls == cl, cl, "_rest_")
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
  colnames (enrich) <- c ('cluster', 'group', 'subgroup', 'fisher.test.pvalue')
  write.csv (enrich, sprintf ("%s-cluster-enrichment.csv", type), row.names=FALSE)
}
