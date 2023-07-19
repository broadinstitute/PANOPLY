args = commandArgs(TRUE)


annot_of_comparison = args[1]

sankey_tar = args[2]
label = args[3]
primary_dataype_label = args[4]

pwd=getwd()
rmarkdown::render("/prot/proteomics/Projects/PGDAC/src/nmf_rmd.rmd",
                  params = list(title = paste0("NMF Report - ", label),
                                sankey_tar = sankey_tar,
                                label = label,
                                annot_of_comparison = annot_of_comparison,
                                primary_dataype_label = primary_dataype_label),
                  output_file = file.path(pwd,paste0(label,"_nmf_rmd.html")))
