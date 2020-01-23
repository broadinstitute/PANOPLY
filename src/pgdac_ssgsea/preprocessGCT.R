#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("optparse"))
suppressPackageStartupMessages(p_load("glue"))

options( warn = -1, stringsAsFactors=F )

# specify command line arguments
option_list <- list(
  make_option( c("-i", "--input"), action='store', type='character',  dest='gct.str', help='Path to input GCT file.'),
  make_option( c("-l", "--level"), action='store', type='character',  dest='level', help="Mode of report, 'ssc' - single-site-centric, 'gc' - gene-centric, 'gcr' - gene-centric-redundant.", default='ssc'),
  make_option( c("-t", "--id_type"), action='store', type='character',  dest='id.type', help="Notation of site-ids: 'sm' - Spectrum Mill; 'wg' - Web Gestalt; 'ph' - Philosopher", default='sm'),
  make_option( c("-o", "--id_type_out"), action='store', type='character',  dest='id.type.out', help="Type of site id for output: 'uniprot', 'refseq', 'seqwin', 'psp' (psp not implemented yet).", default='uniprot'),
  make_option( c("-a", "--acc_type_in"), action='store', type='character',  dest='acc.type', help="Type of accession number in 'rid' object in GCT file (uniprot, refseq, symbol).", default='refseq'),
  make_option( c("-s", "--seqwin_column"), action="store", type='character', dest='seqwin.col', help="Column containing flanking sequences, separated by '|'. Only relevant if '--id_type_out' = 'seqwin'", default='VMsiteFlanks'),
  make_option( c("-g", "--gene_symbol_column"), action='store', type='character',  dest='gene.col', help="Name of column listing gene names; used for gene centric reports.", default='geneSymbol'),
  make_option( c("-v", "--sgt_column"), action='store', type='character',  dest='SGT.col', help="Column used to collpase subgroup-top (SGT) reports.", default='subgroupNum'),
  make_option( c("-d", "--localized"), action='store', type='logical',  dest='loc', help="CAUTION: it is NOT RECOMMENDED to set this flag to FALSE. If TRUE only fully localized sites will be considered." , default=TRUE),
  make_option( c("-m", "--mode"), action='store', type='character',  dest='mode', help="Determines how multiple sites per gene will be combined: sd - most variable (standard deviation) across sample columns; SGT - subgroup top: first subgroup in protein group (Spectrum Mill); abs.max - for log-transformed, signed p-values" , default='median'),
  make_option( c("-r", "--residue"), action='store', type='character',  dest='mod.res', help='Modified residues, e.g. "S|T|Y" or "K".', default='S|T|Y'),
  make_option( c("-p", "--ptm"), action='store', type='character',  dest='mod.type', help='Type of modification, e.g "p" or "ac".', default = 'p'),
  make_option( c("-u", "--preprocess_gct"), action='store', type='logical',  dest='preprocess.gct', help='If FALSE nothing will be done; probably needed for to make this step optional in a FireCLoud WDL.', default = FALSE),
  make_option( c("-z", "--libdir"), action='store', type='character',  dest='libdir', help='Folder to source from', default = 'c:/Users/karsten/Dropbox/Devel/PGDAC/src/pgdac_ssgsea/')
)

# parse command line parameters
opt <- parse_args( OptionParser(option_list=option_list) )

source(glue("{opt$libdir}/pgdac_ptmsea_functions.R"))


###################################################
##       run the function
out <- preprocessGCT(gct.str = opt$gct.str,
              level=opt$level,
              id.type=opt$id.type,
              id.type.out=opt$id.type.out,
              acc.type=opt$acc.type,
              seqwin.col=opt$seqwin.col,
              gene.col=opt$gene.col,
              SGT.col=opt$SGT.col,
              loc=opt$loc,
              mode=opt$mode,
              mod.res = opt$mod.res,
              mod.type=opt$mod.type,
              appenddim=F,
              preprocess.gct=opt$preprocess.gct
)
writeLines(out, con="fn.out")
