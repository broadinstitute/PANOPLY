library(cmapR)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
nmf_tar = args[1]
ome_type = args[2]


untar(nmf_tar) #untar nmf results

#locate K_. directory
subdir=list.dirs(".", recursive=FALSE)
kdir=list.files(subdir, pattern="^K_.$", full.names = TRUE)
# pull out the driver feature gcts for each cluster
ome_gcts_str = list.files(kdir, "*_data_matrix_C._n.+?x.+?.gct", full.names = TRUE)

ome_gcts = c()
for (gct_str in ome_gcts_str) {
  gct = parse_gctx(gct_str)
  if (gct_str==ome_gcts_str[1]) { #if we're on the first gct
    full_gct = gct #save it flatly
  } else { #otherwise
    if(sum(!(gct@rid %in% full_gct@rid))>0) { # if there are non-duplicate rows
      gct_sub = subset_gct(gct, which(!(gct@rid %in% full_gct@rid))) #remove duplicate rows
      full_gct = merge_gct(full_gct, gct_sub) # and merge to the full gct
    }
  }
}

# rename original NMF column
full_gct@cdesc = rename(full_gct@cdesc,
                        NMF.cluster.membership_og = NMF.cluster.membership,
                        NMF.consensus_og = NMF.consensus,
                        NMF.consensus.core_og = NMF.consensus.core,
                        NMF.consensus.mapped_og = NMF.consensus.mapped)

# create gct file
write_gct(full_gct, paste0(ome_type, "-data_matrix_allClusters_n", length(full_gct@cid),'x',length(full_gct@rid), ".gct"))
