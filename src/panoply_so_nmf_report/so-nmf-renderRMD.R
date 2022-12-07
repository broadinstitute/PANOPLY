args = commandArgs(TRUE)

nmf_tar = args[1]
sankey_tar = args[2]
label = args[3]

pwd=getwd()
rmarkdown::render("/prot/proteomics/Projects/PGDAC/src/so_nmf_rmd.rmd",
                  params = list(title = paste0('Single -Ome NMF Comparisons - ', label),
                                nmf_tar = nmf_tar,
                                sankey_tar = sankey_tar,
                                label = label),
                  output_file = file.path(pwd,paste0(label,"_so_nmf_rmd.html")))
