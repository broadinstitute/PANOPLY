
source ('preamble.r')
Source ('misc.r')


## process mutations.tsv table and extract relevant info
## processing is dictated by the format of the mutations.tsv table

d <- read.table (file.path (data.dir, 'mutations.tsv'), 
                 header=TRUE, sep='\t', check.names=FALSE, quote="\"", as.is=TRUE)
d.info <- d [, 1:12]
d <- d [, 13:ncol(d)]
colnames (d) <- sapply (colnames (d), function (x) substr (x, 6, nchar(x)))
d [d==""] <- 0


for (i in 1:ncol(d)) 
  d[,i] <- as.numeric (sapply (d[,i], function (x) ifelse (x==0 || is.na(x), x, 
                                                           length (strsplit (x, split=',')[[1]]))))


accession <- d.info[,'accession_number']
gene <- sapply (d.info[,'gene_name'], function (x) strsplit (x, split=';')[[1]][1])

data <- data.frame (accession.number=accession, gene.name=gene, d, check.names=FALSE)

# remove proteins with missing gene names
data <- data [ !is.na (data[,'gene.name']), ]
# deal with duplicate symbols
repeated.index <- repeated (data[,'gene.name'])
data.unique <- data [ !repeated.index, ]

data.rep <- data [ repeated.index, ]
rep.genes <- data.rep [,'gene.name']
data.rep.chosen <- NULL
for (g in (unique (rep.genes))) {
  data.subset <- data.rep [g==rep.genes,]
  # union of mutations in all rows
  row <- c ( as.character (data.subset[1,1]), as.character (data.subset[1,2]),
             apply (data.subset[,c(-1,-2)], 2, function (x) ifelse (sum(x)>0, 1, 0)) )
  data.rep.chosen <- rbind (data.rep.chosen, row)
}
colnames (data.rep.chosen) <- colnames (data.rep)
data.final <- rbind (data.unique, data.rep.chosen)
# write out results
write.csv (data.final, 'mutated-genes.csv', row.names=FALSE, quote=FALSE)


## binary table
d.binary <- data.final [, c(-1,-2)]
d.binary [ d.binary > 1] <- 1
data.binary <- cbind (data.final[,1:2], d.binary)
# write out results
write.csv (data.binary, 'mutated-genes-binary.csv', row.names=FALSE, quote=FALSE)



## determine top N most frequently mutated genes
for (i in 3:ncol(data.binary)) data.binary [,i] <- as.numeric (data.binary[,i])
mut.n <- 25
mut.count <- apply (data.binary[,c(-1,-2)], 1, sum, na.rm=TRUE)
mut.topN <- data.binary [ mut.count >= mut.n, ]
write.csv (mut.topN, paste ('mutated-genes-in', mut.n, '.csv', sep=''), 
           row.names=FALSE, quote=FALSE)


