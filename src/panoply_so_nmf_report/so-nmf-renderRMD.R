args = commandArgs(TRUE)

so_nmf_tar = args[1]
sankey_tar = args[2]
label = args[3]
mo_nmf_tar = args[4]

pwd=getwd()
rmarkdown::render("/prot/proteomics/Projects/PGDAC/src/so_nmf_rmd.rmd",
                  params = list(title = paste0('Single -Ome NMF Comparisons - ', label),
                                so_nmf_tar = so_nmf_tar,
                                mo_nmf_tar = mo_nmf_tar,
                                sankey_tar = sankey_tar,
                                label = label),
                  output_file = file.path(pwd,paste0(label,"_so_nmf_rmd.html")))
