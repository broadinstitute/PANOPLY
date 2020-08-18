# ### Pipeline parameters
# # defaults shown below
# # can be overridden by redefining the appropriate variable in the
# # parameter redefinition section
# ###
# 
# 
# ## Label type for MS experiment (set.label.type must be called to initialize variables)
# label.type <- 'TMT10'   # alternatives: iTRAQ4, TMT10.126
# set.label.type (label.type) 
# 
# ## QC
# # QC status can be indiated using a separate cls file,
# # or included in the experiment design file with column name qc.col;
# # if neither is found, all samples are marked as qc.pass.label
# # if both exist, the experiment design file supercedes
# sampleQC.cls <- NULL
# qc.col <- 'QC.status'
# qc.pass.label <- 'QC.pass'
# 
# 
# ## Output precision for gct tables
# ndigits <- 5    
# 
# 
# ## Missing values and filtering
# na.max <- 0.7                  # maximum allowed NA values (per protein/site/row), can be fraction or integer number of samples
# min.numratio.fraction <- 0.2   # fraction of samples in which min. numratio should be present to retain protein/phosphosite
# sample.na.max <- 0.8           # maximum allowed fraction of NA values per sample/column; pipeline error if violated
# nmiss.factor <- 0.5            # for some situations, a more stringent condition is needed
# sd.filter.threshold <- 0.5     # SD threshold for SD filtering
# clustering.sd.threshold <- 2   # threshold for filtering data before consensus clustering
# clustering.na.threshold <- 0.5 # max fraction of missing values for clustering; rest are imputed
# apply.SM.filter <- TRUE        # if TRUE, apply numRatio based filter (use TRUE if input is SM ssv)
# 
# 
# ## Normalization
# norm.method <- '2comp'         # options: 2comp (default), median, mean
# alt.method <- 'median'         # alt.method for comparison -- filtered datasets not generated
# if (norm.method == alt.method) alt.method <- NULL
#                                # ignored if alt.method is NULL, or is identical to norm.method
# 
# 
# ## Gene mapping
# # gene mapping not needed -- use SM geneSymbol (but map to official symbols for CNA analysis)
# official.genesyms <- 'gene-symbol-map.csv'
# gene.id.col <- 'geneSymbol'
# protein.gene.map <- 'RefSeq-GeneName-Map-20170701.txt'
# # policy for combining/collapsing duplicate gene names -- uncomment appropriate line to use
# # duplicate.gene.policy <- ifelse (type == 'phosphoproteome', 'median', 'maxvar')  
# duplicate.gene.policy <- 'maxvar'
# 
# ## RNA related
# rna.output.prefix <- 'rna-seq'  # output prefix for tables creates during RNA analysis
# rna.sd.threshold <- 1           # for variation filter (set to NA to disable)
# 
# ## CNA/parallelism related
# pe.max.default <- 250           # default maximum processors/jobs
# 