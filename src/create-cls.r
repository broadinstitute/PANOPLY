

source ('preamble.r')



create.cls <- function (out.prefix, expt.design, type.cls) {
  
  # get proper column names and class information
  tumor.info <- read.csv (expt.design)
  rownames (tumor.info) <- tumor.info[,'Sample.ID']

  sample.order <- scan (paste (master.prefix, 'gct', sep='.'),
                        what='char', nlines=1, skip=2, sep='\t')[c(-1,-2)]
  tumor.group <- tumor.info [ sample.order, 'Subgroup' ]
  
  # global class labels
  write.cls (tumor.group, paste (out.prefix, '-group.cls', sep=''))
  
  # class labels separated by uni/bi-modal groups
  cls <- read.cls (type.cls)
  for (type in c ('unimodal', 'bimodal')) {
    t.cls <- cls == type
    if (any (t.cls)) {
      out <- paste (out.prefix, '-', type, sep='')
      # write all class labels
      write.cls (tumor.group[t.cls], paste (out, '-group.cls', sep=''))
    }
  }  
}



# different "type"s result in the same cls files, 
# but different file prefix is convenient for batch scripts
create.cls (type,
            expt.design=file.path (data.dir, 'exptdesign.csv'),
            type.cls=file.path (data.dir, bimodal.cls))

