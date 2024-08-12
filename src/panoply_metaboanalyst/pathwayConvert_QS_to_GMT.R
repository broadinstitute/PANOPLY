library(ActivePathways) # GMT editing files
library(cmapR)

# read in pathway QS file
pathway = "smpdb_pathway"
df = qs::qread(glue('{pathway}.qs'))

# mapping to HMDB IDs
compound_map = qs::qread('compound_db.qs')
map_to_hmdb=TRUE

gmt = list()
for (i in 1:dim(df)[1]) {
  id = df$id[[i]]
  gene_vec = strsplit(df$member[[i]], '; ')[[1]]
  if (map_to_hmdb) { # if we want to overwrite with HMDB IDs
    gene_vec = compound_map[match(gene_vec, compound_map$name), 'hmdb_id'] # overwrite with HMDB IDs
  }
  if (grepl('^.+?"(http:[^ ]+?)".+?(http:[^ ]+?).+?$', df$reference[[i]])) {
    name = gsub('^.+?"(http:[^ ]+?)".+?(http:[^ ]+?).+?$', '\\2', df$reference[[i]])
  } else {  name = gsub('^.+?"(http:[^ ]+?)".+?$', '\\1', df$reference[[i]]) }
  gmt[[id]] = list(id = id,
                   name = name, # pull reference URL out of reference column
                   genes = gene_vec)
}
# enforce GCT class
class(gmt) <- 'GMT'
# is.GMT(gmt)

# write GCT object
write.GMT(gmt, glue('{pathway}.gmt'))


# file.copy('smpdb_pathway.gmt', '/opt/input/', overwrite=T)



