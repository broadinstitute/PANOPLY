library(glue)

dir_save = '~/Git/panoply-sandbox/src/panoply_metaboanalyst/pathway_db/'


# Pathway Download ####
## pathway database paths are exposed in the MetaboAnalystR GitHub
## the website, however, is more updated than the GitHub


## in general, the path is exposed at
## https://github.com/xia-lab/MetaboAnalystR/blob/master/R/general_lib_utils.R#L8


## Metabolite Set Enrichment Databases ####

## path is exposed at MetaboAnalystR/R/enrich_mset.R
## https://github.com/xia-lab/MetaboAnalystR/blob/master/R/enrich_mset.R#L57
## in the SetCurrentMsetLib() function

## library names are exposed at MetaboAnalystR/man/SetCurrentMsetLib.Rd
## note that 'self' is for user-uploaded data
## https://github.com/xia-lab/MetaboAnalystR/blob/master/man/SetCurrentMsetLib.Rd
# ## unfortunately, this list seems to be out of date with the website
# msets = c(#"self",
#   "kegg_pathway", "smpdb_pathway", "blood", "urine", "csf", "snp", "predicted", "location", "drug")

## pathway names can also be pulled by walking through the web portal for Enrichment Analysis
## https://dev.metaboanalyst.ca/Secure/enrichment/EnrichParamView.xhtml
## and checking the pathway listed in SetCurrentMsetLib() from the R Command History
msets = c(
  'smpdb_pathway',
  'kegg_pathway', # 	81 metabolite sets based on KEGG human metabolic pathways (Nov. 2025)
  'kegg_pathway_2023', # 80 metabolite sets based on KEGG human metabolic pathways (Dec. 2023) [Deprecated]
  'drug', # 461 metabolite sets based on drug pathways from SMPDB.
  'RaMP_pathway', #	3425 metabolite sets based on RaMP-DB pathways.
  'lipid_pathway', # 	817 entries integrating PathBank, Reactome, WikiPathways, and KEGG (Nov. 2025).
  'kegg_gutbac', # 	123 entries based on KEGG metabolic pathways of common human gut bacteria 
  'kegg_gutbachsa' # 139 KEGG pathways integrating human (host) and gut bacterial metabolism.
)

for (mset in msets) {
  fn = glue("https://www.metaboanalyst.ca/resources/libs/msets/{mset}.qs")
  download.file(fn, file.path(dir_save, basename(fn)))
  # df = qs::qread(file.path(pathway.dir,glue('{pathway}.qs')))
}

## Excluded Pathways ####
# Disease signatures
# 
# Blood
# 480 metabolite sets reported in human blood.
# 
# Urine	385 metabolite sets reported in human urine.
# 
# CSF	174 metabolite sets reported in human cerebral spinal fluid (CSF).
# 
# Feces	67 metabolite sets reported in human feces.
# Chemcial structures
# 
# Super-class
# 39 super chemical class metabolite sets or lipid sets
# 
# Main-class	617 main chemical class metabolite sets or lipid sets
# 
# Sub-class	1250 sub chemical class metabolite sets or lipid sets
# Other types
# 
# SNPs
# 4,598 metabolite sets based on their associations with SNPs loci.
# 
# Exposure	62 metabolite sets based on dietary and chemical exposures.
# 
# Predicted	912 metabolic sets predicted to change in the case of dysfunctional enzymes.
# 
# Locations	78 metabolite and lipid sets based on organ, tissue, and subcellular localizations.




## Integrated Analysis Pathway Databases ####

## path is exposed at R/enrich_integ.R
## https://github.com/xia-lab/MetaboAnalystR/blob/master/R/enrich_integ.R#L141
## we are interested in the 'integ' or integrated set, containing both metabolites and genes
## we are only interested in the human database ('hsa.qs')

for (org in c('hsa')) {
  fn = glue("https://www.metaboanalyst.ca/resources/libs/kegg/jointpa/integ/{org}.qs")
  download.file(fn, file.path(dir_save, basename(fn)))
}

# Additional Databases ####
## the compound_db.qs object contains mappings between compound names and various other IDs
## the filename is exposed at
## https://github.com/xia-lab/MetaboAnalystR/blob/master/R/util_approx.R
## there appears to be a lipid-specific database, as well as a 'master' database


fn = 'https://www.metaboanalyst.ca/resources/libs/master_compound_db'
download.file(fn)


# Testing Pathway DBs ####
for (mset in msets) {
  df = qs::qread(file.path(dir_save,glue('{mset}.qs')))
  # all have the same general structure
  print( names(df) )
  # all use Compound Names as their identifiers
  print( df$member[1] )
} 

