#!/usr/bin/env Rscript
options( warn = -1 )
args <- commandArgs(trailingOnly=T)

## get arguments
mut.vcf <- args[1]
junc.bed <- args[2]
##fus.bed <- args[3]
id <- args[3]
annotation.zip <- args[4]

## used to extract annotation files
tmp.dir <- '/tmp/' 
##logfile <- args[7]

## ######################################################################################################
##
##                run customProDB on a single sample
##
## mut.vcf           character, mutation calls in VCF format
## junc.bed          NULL or character, junctions in BED format
## fus.bed           NULL or character, jusions in BED format (not supported yet)
## id                character, sample id. Will be used as filename prefic for
##                   result files
## annotation.zip    Zip-archive containing annotation files downloaded as described
##                   in the customProDB R-Vignette
## tmp.dir           character, path to a writeable folder. Used to extract annotation files
##
## ######################################################################################################
run_cpdb <- function(mut.vcf, junc.bed=NULL, fus.bed=NULL, id='test', annotation.zip='./',  tmp.dir='.'){

    ## packages
    require(customProDBBI)
    require(WriteXLS)
    require(BSgenome.Hsapiens.UCSC.hg19)

    ## prepare log file
    logfile=paste('cpdb_', id, '.log', sep='')

    start.time <- Sys.time()
    cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'run_cpdb\'--\n\n', file=logfile)
    cat('## parameters\nmutation file:', mut.vcf, '\njunction file:', junc.bed, '\nfusion file:', fus.bed, '\nannotation file:', annotation.zip, '\ntemp-dir:', tmp.dir, '\nlog file:', logfile, '\n', file=logfile, append=T)


    ## ###########################################
    ## temporary id
    fn.tmp <- paste(paste(letters[sample(26, 5)], collapse=''), paste(sample(100,5), collapse=''), sep='')

    ## ###########################################
    ## unzip annotation file
    ## prepare tmp folder
    tmp.dir <- paste(tmp.dir, fn.tmp, sep=ifelse( substr(tmp.dir, nchar(tmp.dir),nchar(tmp.dir)) == '/' , '', '/'))
    dir.create( tmp.dir  )

    ## unzip
    cat('\n## Extracting annotation files to', tmp.dir, '\n', file=logfile, append=T)
    unzip(annotation.zip, exdir=tmp.dir)


    ## ###########################################
    ## filename for output files
    fn.tmp <- id

    ## ###########################################
    ## import annotation files
    cat('\n## Importing annotation files\n', file=logfile, append=T)
    load(paste(tmp.dir, "exon_anno.RData", sep="/"))
    load(paste(tmp.dir, "ids.RData", sep="/"))
    load(paste(tmp.dir, "proseq.RData", sep="/"))
    load(paste(tmp.dir, "procodingseq.RData", sep="/"))
    load(paste(tmp.dir, "splicemax.RData", sep="/"))
    txdb <- loadDb( paste(tmp.dir, "txdb.sqlite", sep='/') )


    ## ########################################################
    ## check whether chromosome starts with 'chr'
    ## fix, if it doesn't...
    if(length(grep('^chr', read.delim(mut.vcf, nrow=1, comment.char='#', header=F)[1])) == 0 ){
        vcf <- readLines(mut.vcf)
        vcf.head.idx <- grep('^#', vcf)
        vcf.head <- vcf[vcf.head.idx]
        vcf <- paste('chr', vcf[-vcf.head.idx], sep='')
        mut.vcf <- sub('\\.vcf', '_fix.vcf', mut.vcf)
        writeLines( c(vcf.head, vcf), con=mut.vcf)
    }

    ## prepare output directory

    ## #########################################################
    ##
    ##                 splice junction database
    ##
    ## ##########################################################
    if(!is.null(junc.bed)){

        cat('\n## Junction pipeline\n', file=logfile, append=T)

        ## import junctions
        cat('Importing splice junctions\n', file=logfile, append=T)
        jun= Bed2Range(junc.bed)

        ## determine junction type
        cat('Annotating splice junctions\n', file=logfile, append=T)
        jun_type <- JunctionType(jun, splicemax, txdb, ids)

        ## extract 'regular' chromosomes
        jun_type <- jun_type[grep('chr([0-9]{1,2}|X|Y)', jun_type$seqnames, value=F),]

        ## write fasta
        cat('Exporting junction fasta\n', file=logfile, append=T)
        OutputNovelJun <- OutputNovelJun(jun_type, Hsapiens, paste(fn.tmp, '_junctions.fasta', sep=''),  proteinseq)
    }

    ## ##########################################################
    ##
    ##                   gene fusion database
    ## ToDo
    ## ###########################################################





    ## ###########################################################
    ##
    ##                 mutation database
    ##
    ## ###########################################################
    if(!is.null(mut.vcf)){
        cat('\n## VCF pipeline\n', file=logfile, append=T)


        ## ######################
        ## import variants
        cat('Importing VCF\n', file=logfile, append=T)
        vcf <- InputVcf(mut.vcf)

        ## separate INDELS from SNVs
        cat('Extracting indels... ', file=logfile, append=T)
        indel.idx <- which(nchar(values(vcf[[1]])$ALT )  != nchar(values(vcf[[1]])$REF))
        vcf.indel <- vcf[[1]][ indel.idx ]
        cat(length(vcf.indel),'indels found\n' ,file=logfile, append=T)

        cat('Extracting SNVs...', file=logfile, append=T)
        snv.idx <- which(nchar(values(vcf[[1]])$ALT)  == nchar(values(vcf[[1]])$REF))
        vcf.snv <- vcf[[1]][ snv.idx ]
        cat(length(vcf.snv),'SNVs found\n' ,file=logfile, append=T)

        ## #######################################
        ##                SNV database
        ## #######################################
        if(length(vcf.snv) > 0){

            cat('\n## parsing SNVs\n' ,file=logfile, append=T)

            ## ##########################################
            ##              predict locations
            ## ##########################################
            cat('predicting SNV locations\n' ,file=logfile, append=T)
            SNVloc <- Varlocation(vcf.snv, txdb, ids)

            ## distribution of locations
            tab.snv <- table( SNVloc[, "location"] )

            ## ###########################################
            ##         extract coding variants
            ## ###########################################
            ## position of snvs in coding sequence

            cat('Extracting SNVs in coding regions\n' ,file=logfile, append=T)
            postable.snv <- Positionincoding( vcf.snv, exon  )

            ## ###########################################
            ##        remove DNPs for the time being
            ## ############################################
            dnp.idx <- which(nchar(postable.snv$refbase) > 1)
            ##View(postable.snv)
            postable.snv.org <- postable.snv

            cat('Removing DNVs\n' ,file=logfile, append=T)
            if(length(dnp.idx) > 0)
                postable.snv <- postable.snv[-dnp.idx, ]

            ## ############################################
            ##              predict nsSNV
            ## ############################################
            txlist <- unique(postable.snv[, "txid"])

            cat('Extracting protein coding sequences affected by SNVs\n' ,file=logfile, append=T)
            ## extract protein coding sequences affected by SNVs
            codingseq <- procodingseq[ procodingseq[, "tx_id"] %in% txlist, ]
            
            cat('Predicting nsSNVs\n' ,file=logfile, append=T)
           
            ## predict non-synonymous SNVs
            mtab <- aaVariation(postable.snv, codingseq)

            ## ##########################
            ## export variant database
            ## ##########################
            cat('Exporting nsSNVs fasta\n' ,file=logfile, append=T)
			if(nrow(mtab) > 0)
            	OutputVarproseq(mtab, proteinseq, paste(fn.tmp, '_nsSNV.fasta', sep=''), ids, RPKM=NULL)
			else
				cat('', file=paste(fn.tmp, '_nsSNV.fasta', sep=''))

            ## ##########################
            ## export table
            codingseq=data.frame(codingseq)
            postable.snv <- data.frame(postable.snv.org)
            mtab=data.frame(mtab)
            tab.snv=data.frame(tab.snv)

            ##WriteXLS(c('postable.snv', 'mtab', 'tab.snv'), ExcelFileName=paste(fn.tmp, '_nsSNV_INFO.xlsx', sep=''), SheetNames=c('SNV coding', 'non-synon', 'location table'), FreezeRow=1, AutoFilter=T,BoldHeaderRow=T )
        }
        ## #######################################################
        ##                    INDEL database
        ## ########################################################
        if(length(vcf.indel)){

            ## #################################
            ##              predict locations
            ## ##################################
            cat('\n## parsing indels\n' ,file=logfile, append=T)

            cat('Extracting indels in coding regions\n' ,file=logfile, append=T)
            postable.indel <- Positionincoding( vcf.indel, exon  )
            txlist.indel <- unique(postable.indel[, "txid"])

            cat('Predicting indel locations\n' ,file=logfile, append=T)
            indelloc <- Varlocation(vcf.indel, txdb, ids)

            tab.indel <- table( indelloc[, "location"] )

            ## extract protein coding sequences affected by SNVs
            codingseq.indel <- procodingseq[ procodingseq[, "tx_id"] %in% txlist.indel, ]

            ## export
            cat('Exporting indel fasta\n' ,file=logfile, append=T)
			if(nrow(postable.indel) > 0)
            	Outputaberrant(postable.indel, coding=codingseq.indel, proteinseq=proteinseq, outfile=paste(fn.tmp, '_indel.fasta', sep=''), ids, RPKM=NULL)
			else 
				cat('', file=paste(fn.tmp, '_nsSNV.fasta', sep=''))


            ## ##########################
            ## export table
            postable.indel <- data.frame(postable.indel)
            codingseq.indel <- data.frame(codingseq.indel)
            indelloc <- data.frame(indelloc)
            tab.indel <- data.frame(tab.indel)

            ## writing out the nucleotide sequence (codingseq.indel) causes perl to throw an error due to Excel's character limit per cell
            ##WriteXLS(c('postable.indel', 'indelloc', 'tab.indel'), ExcelFileName=paste(fn.tmp, '_indel_INFO.xlsx', sep=''), SheetNames=c('indel coding', 'indel location', 'location table'), FreezeRow=1, AutoFilter=T,BoldHeaderRow=T )
        }
    }

    ## ############################################
    ## delete tmp dir
    ##saveDb(txdb,  paste(tmp.dir, "txdb.sqlite", sep='/') )

    ##cat('\n## Cleaning up...' ,file=logfile, append=T)
    ##file.remove(dir(tmp.dir, full.names=T))
    ##unlink(tmp.dir)

    cat('done\n' ,file=logfile, append=T)
    cat('\nTotal time:', paste0(format(Sys.time() - start.time)), file=logfile, append=T)

    return(0)
} ## end function 'run_cpdb'


## #####################################################
##                          run
##res <- run_cpdb(mut.vcf=mut.vcf, junc.bed=junc.bed, fus.bed=fus.bed, id=id, annotation.zip=annotation.zip, tmp.dir=tmp.dir )
res <- run_cpdb(mut.vcf=mut.vcf, junc.bed=junc.bed, id=id, annotation.zip=annotation.zip, tmp.dir=tmp.dir )



