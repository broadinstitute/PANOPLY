
### Pipeline parameters -- uncomment as needed (or generate from code)


# ## Label type for MS experiment (set.label.type must be called to initialize variables)
# label.type <- 'TMT10'   # alternatives: iTRAQ4
# set.label.type (label.type) 
# 
# ## Sample replicate indicator
# # Sample.IDs MUST be unique in the expt.design.file; duplicate samples should have the same 
# # sample names, but include this replicate.indicator, followed by a unique suffix: <name>REP1)
# replicate.indicator <- '.REP'
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
# # nmiss.plex <- 0.25
# # n.plex <- length (unique (read.csv (expt.design.file)[,'Experiment']))
# # na.max <- ceiling (n.plex * length(plex.channels) * nmiss.plex)                  
# #                               # must be present in at least nmiss.plex fraction of the experiments (plexes)
# na.max <- 0.7                 # maximum allowed NA values, can be fraction or integer number of samples
# nmiss.factor <- 0.5           # for some situations, a more stringent condition is needed
# sd.filter.threshold <- 0.5    # SD threshold for SD filtering
# apply.SM.filter <- TRUE       # if TRUE, apply numRatio based filter (use TRUE if input is SM ssv)
# 
# 
# ## Normalization
# norm.method <- 'median'        # options: 2comp (default), median, mean
# alt.method <- 'median'         # alt.method for comparison -- filtered datasets not generated
# if (norm.method == alt.method) alt.method <- NULL
# # ignored if alt.method is NULL, or is identical to norm.method
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
# ## Project
# # data source -- for managing some operations (esp related to sample IDs and names)
# #  [all current options listed below -- uncomment only one]
# # use project.name to manage project specific processing and options
# project.name <- 'default'
# # project.name <- 'cptac2.tcga'
# 
# ## Disease
# # disease setting is used to set disease specific options and
# # run appropriate CNA subsets by creating run-cna-analysis-<disease>.r;
# # uncomment one (or none) below
# # disease <- 'MEDULLO'
# # disease <- 'BRCA'
