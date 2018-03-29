
source ('config.r')

best.clus <- read.table ( sprintf ("%s.bestclus.txt", type), header=TRUE)
best.clus <- best.clus [,-3]
colnames (best.clus) <- c ('Sample.ID', 'Cluster')

write.csv (best.clus, sprintf ("%s-bestclus.csv", type), row.names=FALSE, quote=FALSE)

# add groups file to config.r for running association analysis
write ( sprintf ("assoc.subgroups <- '%s-bestclus.csv'", type), file="config.r", append=TRUE )
