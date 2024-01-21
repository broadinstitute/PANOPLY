args = commandArgs(TRUE)


rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  #### NMF Outputs ####
  make_option( c("-n", "--nmf_results"), action='store', type='character',  dest='nmf_results', help='Rdata object with NMF Clustering results (res.rank object).'),
  make_option( c("-r", "--rank_top"), action='store', type='numeric',  dest='rank_top', help='Best number of clustering / rank.'),
  make_option( c("-e", "--expr_comb"), action='store', type='character',  dest='expr_comb', help='GCT file with combined expression data.'),
  make_option( c("-f", "--expr_comb_nn"), action='store', type='character',  dest='expr_comb_nn', help='GCT file with combined (non-negative) expression data, used for NMF.'),
  make_option( c("-p", "--nmf_parameters"), action='store', type='character',  dest='nmf_parameters', help='Rdata file with parameters (opt object) from panoply_nmf.'),
  #### Postprocessing Outputs ####
  make_option( c("-t", "--postprocess_tar"), action='store', type='character',  dest='postprocess_tar', help='Tar file containing figures and analyses for NMF results.'),
  make_option( c("-q", "--postprocess_parameters"), action='store', type='character',  dest='postprocess_parameters', help='Rdata file with parameters (opt object) from panoply_nmf_postprocess.'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='label', help='Label associated with this run.'),  # default = 2),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)

opt <- parse_args( OptionParser(option_list=option_list) )

# # testing locally
# opt=list()
# opt$nmf_results="/Users/wcorinne/Downloads/panoply_nmf_report/inputs/nmf_res.Rdata"
# opt$nmf_parameters="/Users/wcorinne/Downloads/panoply_nmf_report/inputs/nmf_opt.Rdata"
# opt$postprocess_parameters="/Users/wcorinne/Downloads/panoply_nmf_report/inputs/nmf_postprocess_opt.Rdata"
# opt$label="luad_test"
# opt$rank_top=6
# opt$expr_comb="/Users/wcorinne/Downloads/panoply_nmf_report/inputs/luad_test_combined_n163x51372.gct"
# opt$expr_comb_nn="/Users/wcorinne/Downloads/panoply_nmf_report/inputs/luad_test_combinedNonNegative_n163x102744.gct"
# opt$postprocess_tar="/Users/wcorinne/Downloads/panoply_nmf_report/inputs/NMF_results.tar.gz"
# opt$lib_dir="~/Git/panoply-sandbox/src/panoply_nmf_report"

# print parameters
cat("\n\n####################\nNMF POST-PROCESSING PARAMETERS\n\n")
print(opt)
save(file = "nmf_postprocess_opt.Rdata", opt)
cat("####################\n\n")


pwd=getwd()
rmarkdown::render(file.path(opt$lib_dir,"nmf_rmd.rmd"),
                  params = list(title = paste0("NMF Report - ", opt$label),
                                nmf_results = opt$nmf_results,
                                nmf_parameters = opt$nmf_parameters,
                                postprocess_parameters = opt$postprocess_parameters,
                                label = opt$label,
                                gct_expr_comb_file = opt$expr_comb,
                                gct_expr_comb_nn_file = opt$expr_comb_nn,
                                rank_top = opt$rank_top,
                                postprocess_tar = opt$postprocess_tar),
                  output_file = file.path(pwd,paste0(opt$label,"_nmf_report.html")))


