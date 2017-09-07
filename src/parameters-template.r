


### Pipeline parameters -- uncomment as needed (or generate from code)

# 
#
# ## Label type for MS experiment (set.label.type must be called to initialize variables)
# label.type <- 'TMT10'   # alternatives: iTRAQ4
# set.label.type (label.type)
#
#
# ## QC
# bimodal.cls <- NULL
# qc.col <- 'QC.status'
# 
# 
# ## Output precision for gct tables
# ndigits <- 5    
# 
# 
# ## Missing values
# nmiss.plex <- 0.25
# n.plex <- length (unique (read.csv (expt.design.file)[,'Experiment']))
# na.max <- ceiling (n.plex * nmiss.plex)                  
# # must be present in at least nmiss.plex fraction of the experiments (plexes)
# nmiss.factor <- 0.5           # for some situations, a more stringent condition is needed
# sd.filter.threshold <- 0.5    # SD threshold for SD filtering
# 
# 
# ## Normalization
# norm.method <- 'median'        # options: 2comp (default), median, mean
# alt.method <- '2comp'          # alt.method for comparison -- filtered datasets not generated
# if (norm.method == alt.method) alt.method <- NULL
# # ignored if alt.method is NULL, or is identical to norm.method
# 
# 
# ## Gene mapping
# # gene mapping not needed -- use SM geneSymbol (but map to official symbols for CNA analysis)
# official.genesyms <- 'gene-symbol-map.csv'
# gene.id.col <- 'geneSymbol'
# # policy for combining/collapsing duplicate gene names -- uncomment appropriate line to use
# # duplicate.gene.policy <- ifelse (type == 'phosphoproteome', 'median', 'maxvar')  
# duplicate.gene.policy <- 'maxvar'
# 
# 
# ## Parallelization
# # mutation and correlation analysis -- max number of parallel jobs
# LSF.mut.jid.max <- 400
# 
# 
# ## Project
# # data source -- for managing some operations (esp related to sample IDs and names)
# #  [all current options listed below -- uncomment only one]
# # use project.name to manage project specific processing and options
# project.name <- 'default'
# # project.name <- 'cptac2.tcga'
# 
