
source ('preamble.r')
Source ('markersel-classify.r')
Source ('stats.r')


class.vectors <- c ('Subgroup', 'cSubgroup')

gct.file <- file.path (norm.dir, paste (master.prefix, '.gct', sep=''))
cls.files <- file.path (norm.dir, paste (type, subset.str, '-', class.vectors, '.cls', sep=''))


# marker selection and classification for each file in cls.files
for (i in 1:length(class.vectors)) {   
  marker.selection.and.classification (gct.file, cls.files[i], paste (class.vectors[i], '-analysis', sep=''), 
                                       id.col=id.col, desc.col=desc.col, gsea=TRUE,
                                       duplicate.gene.policy=duplicate.gene.policy)
}

