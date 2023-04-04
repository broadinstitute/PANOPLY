#library(rtracklayer)
#a <- readGFF("~/Downloads/Homo_sapiens.GRCh38.99.gff3")
#gene_info <- a %>% as.data.frame() %>% dplyr::select(seqid,ID,type, gene_id, Name) %>% filter(type=="gene") %>% dplyr::select(seqid,Name) 

library(dplyr)
library(readr)

# 04/07/2020
# http://useast.ensembl.org/index.html
a <- read.delim("~/Downloads/mart_export.txt",stringsAsFactors = FALSE)
names(a) <- c("gene_id","chromosome_name","hgnc_symbol")
a %>% dplyr::select(hgnc_symbol,chromosome_name) %>% write_tsv("genes.tsv")
