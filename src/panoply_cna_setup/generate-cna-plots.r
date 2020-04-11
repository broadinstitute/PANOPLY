
### Original code from Vanderbilt University. Author: Jing Wang 04212014
## cna_mrna is a matrix only containing 1, 0 and -1. 1 and -1 represent significant positive and negative spearman correlation between mRNA (each row) and CNA (each column)
## cna_protein is a matrix only containing 1, 0 and -1. 1 and -1 represent significant positive and negative spearman correlation between protein (each row) and CNA (each column).
## cna_protein should have the same dimension with cna_mrna. And two matrix should have the row and column names.
## genelocate contains the chromsome name, gene start (bp) and gene end (bp) for each gene
## chromLength contains the chromosome length for each of 24 chromosomes


Plot_cis_trans_effect <- function(cna_mrna,cna_protein,genelocate,chromLength,outputfile){
    
    chrome <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
    genelocate <- genelocate[genelocate[,2] %in% chrome,]
    
    allgene_protein_mrna_cnv <- union(rownames(cna_mrna),colnames(cna_mrna))
    
    allgene_locate <- genelocate[genelocate[,1] %in% allgene_protein_mrna_cnv,]
    
    x <- c(0,chromLength[1:(nrow(chromLength)-1),2])
    chromLength[,3] <- cumsum(as.numeric(x))
    chromLength[,4] <- cumsum(as.numeric(chromLength[,2]))
    
    allChromlen <- chromLength[nrow(chromLength),4]
    
    allgene_locate <- cbind(allgene_locate,0,0)
    
    colnames(allgene_locate)[5:6] <- c("finalstart","finalend")
    
    cat("Calculate gene location in chromosome...\n")
    
    for(i in c(1:nrow(allgene_locate))){
        chr <- allgene_locate[i,2]
        s <- allgene_locate[i,3]
        e <- allgene_locate[i,4]
        cs <- chromLength[chromLength[,1]==chr,3]
        allgene_locate[i,5] <- s+cs
        allgene_locate[i,6] <- e+cs
    }


    ##########Plot Figure############
    
    cat("Plot figure...\n")
    
    #######Identify the common pairs for cna_mrna and cna_protein
    common <- intersect (rownames (cna_mrna), rownames (cna_protein))
    x <- abs(cna_mrna[common,])+abs(cna_protein[common,])
    x[x!=2] <- 0
    ov <- apply(x,2,sum)
    ov <- (-ov)/2

    ######Identify specific pairs for cna_mrna or cna_protein
#     x <- abs(cna_mrna)+abs(cna_protein)
#     x[x==2] <- 0
#     y <- abs(cna_mrna)*x
    y <- abs(cna_mrna)
    spe_mrna <- apply(y,2,sum)
    
#    y <- abs(cna_protein)*x
    y <- abs(cna_protein)
    spe_protein <- apply(y,2,sum)

    maxM <- max(spe_mrna)
    maxP <- max(spe_protein)
    maxO <- min(ov)
    maxS <- max(maxM,maxP)

    png(outputfile,height=480*8,width=480*11,res=300,type='cairo')

    layout(matrix(c(1,2,3,4),2,2),heights=c(2,1))

    par(mar=c(0,4,0,0))

    p <- which(cna_mrna!=0)
    rownum <- nrow(cna_mrna)
    allcnagene <- colnames(cna_mrna)
    allovgene <- rownames(cna_mrna)

    la <- 1
    for(i in c(1:length(p))){
        po <- p[i]
        rowi <- po %% rownum
        if(rowi == 0){
            rowi <- rownum
        }
        coli <- ceiling(po/rownum)
        cnag <- allcnagene[coli]
        ovg <- allovgene[rowi]
        cnagp <- allgene_locate[allgene_locate[,1]==cnag,5]
        ovgp <- allgene_locate[allgene_locate[,1]==ovg,5]
    
        if(length(cnagp)==0 || length(ovgp)==0){
            next
        }
    
        cov <- cna_mrna[rowi,coli]
        color <- ifelse(cov>0,"red","green")
       	if(la==1){
		plot(cnagp,ovgp,xlim=c(0,allChromlen),ylim=c(0,allChromlen),xaxt="n",yaxt="n",frame.plot=F,xlab="",ylab="",pch=20,col=color,cex=0.2)
            axis(side=2,at=(chromLength[,4]-chromLength[,2]/2),labels=chrome)
            abline(h=c(0,chromLength[,4]),v=c(0,chromLength[,4]),col="gray",lty=3)
            la <- la+1
        }else{
            for(u in c(1:length(cnagp))){
                for(v in c(1:length(ovgp))){
                    points(cnagp[u],ovgp[v],pch=20,col=color,cex=0.2)
                }
            }
        }
	}

    par(mar=c(4,4,0,0))

    plot(0,0,xlim=c(0,allChromlen),ylim=c(maxO,maxS),type="n",xaxt="n",frame.plot=F,xlab="",ylab="")
    axis(side=1,at=(chromLength[,4]-chromLength[,2]/2),labels=chrome)

    abline(v=c(0,chromLength[,4]),col="gray",lty=3)

    for(i in c(1:length(spe_mrna))){
        gg <- names(spe_mrna)[i]
        gg_p <- allgene_locate[allgene_locate[,1]==gg,5]
        if(length(gg_p)==0){
            next
        }
        up <- spe_mrna[gg]
        do <- ov[gg]
        for(j in c(1:length(gg_p))){
            points(gg_p[j],up,cex=0.2,type="h",col="blue")
            points(gg_p[j],do,cex=0.2,type="h",col="black")
        }
    }


    ##############Protein#####################
    
    par(mar=c(0,1,0,0))

    p <- which(cna_protein!=0)

    la <- 1
    for(i in c(1:length(p))){
        po <- p[i]
        rowi <- po %% rownum
        if(rowi == 0){
            rowi <- rownum
        }
        coli <- ceiling(po/rownum)
        cnag <- allcnagene[coli]
        ovg <- allovgene[rowi]
        cnagp <- allgene_locate[allgene_locate[,1]==cnag,5]
        ovgp <- allgene_locate[allgene_locate[,1]==ovg,5]
    
        if(length(cnagp)==0 || length(ovgp)==0){
            next
        }
    
        cov <- cna_protein[rowi,coli]
        color <- ifelse(cov>0,"red","green")
        if(la==1){
		plot(cnagp,ovgp,xlim=c(0,allChromlen),ylim=c(0,allChromlen),xaxt="n",yaxt="n",frame.plot=F,xlab="",ylab="",pch=20,col=color,cex=0.2)
            abline(h=c(0,chromLength[,4]),v=c(0,chromLength[,4]),col="gray",lty=3)
            la <- la+1
        }else{
            for(u in c(1:length(cnagp))){
                for(v in c(1:length(ovgp))){
                    points(cnagp[u],ovgp[v],pch=20,col=color,cex=0.2)
                }
            }
        }
    }

    par(mar=c(4,1,0,0))


    plot(0,0,xlim=c(0,allChromlen),ylim=c(maxO,maxS),type="n",xaxt="n",yaxt="n",frame.plot=F,xlab="",ylab="")
    axis(side=1,at=(chromLength[,4]-chromLength[,2]/2),labels=chrome)
    abline(v=c(0,chromLength[,4]),col="gray",lty=3)

    for(i in c(1:length(spe_protein))){
        gg <- names(spe_protein)[i]
        gg_p <- allgene_locate[allgene_locate[,1]==gg,5]
        if(length(gg_p)==0){
            next
        }
        up <- spe_protein[gg]
        do <- ov[gg]
        for(j in c(1:length(gg_p))){
            points(gg_p[j],up,cex=0.2,type="h",col="blue")
            points(gg_p[j],do,cex=0.2,type="h",col="black")
        }
	}


    dev.off()

}
