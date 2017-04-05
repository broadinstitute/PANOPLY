

source ('preamble.r')


sd.threshold <- 1
out.prefix <- 'rna-seq'



d <- read.delim ( file.path (data.dir, 'rnaseq-data.txt'), check.names=FALSE )
data <- d[, -1]
id <- sapply (d[,1], function (x) s <- strsplit (toString (x), split='\\|')[[1]][1])
samples <- colnames (data)

tumor.info <- read.csv ( file.path (data.dir, 'exptdesign.csv'), skip=2, header=F )
samples.order <- as.vector (t (tumor.info[,5:7]) )
run.order <- sapply (samples.order,
                     function (x) {
                       pos <- which (samples %in% x)
                       # missing samples point to NA column add below
                       return ( ifelse (length(pos)==0, ncol(data)+1, pos) )
                     })

data <- cbind (data, na.col=rep(NA, nrow(data)))  # add NA column at the end -- use for missing samples
data <- cbind (Name=id, Description=id, data[,run.order])
colnames (data) <- c ("Name", "Description", samples.order)

write.gct (data, paste (out.prefix, '.gct', sep=''))


## get PAM 50 genes only
g <- read.table ( file.path (data.dir, 'pam50_annotation.txt'), header=T ) [,1]
pam50.genes <- unlist (lapply (g, toString))
data.pam50 <- data [ data[,'Name'] %in% pam50.genes, ]
write.gct (data.pam50, paste (out.prefix, '-pam50.gct', sep=''))



## filters
# exclude rows with SD less than sd.threshold
sds <- apply (data [, 3:ncol(data)], 1, sd, na.rm=TRUE)
data.subset <- data [ sds > sd.threshold, ]
write.gct (data.subset, paste (out.prefix, '-sdfilter.gct', sep=''))
