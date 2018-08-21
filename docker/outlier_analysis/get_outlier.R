#######################################################################################
##R code for extracting outliers for ovarian cancer data
##input: 
#a).phosphorylation,proteome,copy number and transcriptome expression txt
#b).Group information about subjects' cancer subtypes
##output: 
#a).outlier indicator (for phosphorylation, it's sum of outlier) for each gene, .csv
#b).Outlier indicator for genes within each subtype, .csv
#c).Heatmap of outlier for 4 datasets,arranged by cancer subtypes, .png or .pdf
##Hua Zhou 08/21/2017
#######################################################################################
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

readdata=function(suffix) {
  files=read.table(paste('/input/data/sample_input_', suffix,'.txt',sep='',collapse=''),header=F,sep='\t',comment.char="",stringsAsFactors = F)
  data_raw=list()
  cm_subjs=c()
  cm_genes=c() 
  for ( i in 1:nrow(files)) {
    temp=read.table(paste('/input/data/data/', files[i,1],collapse='',sep=''),header=T,sep='\t',stringsAsFactors = F)
    #check null data counts
    cnt_nan=apply(!(apply(temp[,-1],2,is.na)),1,sum)
    ##remove gene name= null
    temp=temp[temp$gene!='' & cnt_nan>0,]
    temp$gene=unlist(sapply(strsplit(temp$gene,'-'),head,1))
    cm_subjs=union(cm_subjs,names(temp))
    cm_genes=union(cm_genes,temp$gene)
    ##preprocessing data
    ##genename simplified, then if phospho, take first of group, if others,take max
    if (files[i,3]=='phospho')
      temp %>% group_by(gene) %>% summarize_all(funs(max))->temp
    else
      temp %>% group_by(gene) %>% summarize_all(funs(first))-> temp
    if (files[i,6]=='True') {
      med=apply(temp[,-1],1,median,na.rm=T)  
      iqr=apply(temp[,-1],1,IQR,na.rm=T)
      temp[,-1]=(temp[,-1]-med)/iqr } ##get IQR fold change
    data_raw=c(data_raw,list(temp))}
  cm_subjs=c('gene',sort(setdiff(cm_subjs,'gene')))
  ##reformat dataframe by common names (column names)
  for ( i in 1:nrow(files)) {
    cur_data=data_raw[[i]]
    temp=data.frame(matrix(NaN,length(cm_genes),length(cm_subjs)))
    names(temp)=cm_subjs
    temp$gene=cm_genes
    ind1=match(names(cur_data),cm_subjs)
    ind2=match(cur_data$gene,cm_genes)
    temp[ind2[!is.na(ind2)],ind1[!is.na(ind1)]]=cur_data
    data_raw[[i]]=temp}
  return(data_raw) }

## function to get outlier
##median +1.5 IQR, general, for copy number use 1 as cutoff
get_outlier= function(suffix, direction='up') {
  files=read.table(paste('/input/data/sample_input_', suffix,'.txt',sep='',collapse=''),header=F,sep='\t',comment.char="",stringsAsFactors = F)
  data_raw=readdata(suffix)
  for (i in 1:nrow(files)) {
    data=data_raw[[i]]
    cutoff=files[i,5]
    vec=data[,-1]
    switch(direction, 
           up={vec=(vec>=cutoff)},
           down={vec=(vec<=(-cutoff))},
           both={vec=(abs(vec)>=cutoff)})
    
    vec[is.na(vec)]=0
    result=data
    result[,-1]=vec
    res_dir=paste('/output/res/result/',suffix,sep='')
    if (!dir.exists(res_dir))
      dir.create(res_dir,recursive=T)
    filename=paste('/output/res/result/',suffix,'/',files[i,2],'_outlier_', direction,'.txt',sep='')
    if (!file.exists(filename))
      write.table(result,filename,sep='\t',row.names=F,col.names=T)
  } }


plot_outlier=function(suffix,direction,gene_of_interest) {
  files=read.table(paste('/input/data/sample_input_', suffix,'.txt',sep='',collapse=''),header=F,sep='\t',comment.char="",stringsAsFactors = F)
  img_dir=paste('/output/res/result/',suffix,'/img/',direction,sep='')
  if (!dir.exists(img_dir))
    dir.create(img_dir,recursive=T)
  out_mat=list()
  for (i in 1:nrow(files)) {
    out_mat[[i]]=read.table(paste('/output/res/result/', suffix,'/', files[i,2],'_outlier_',direction, '.txt',sep='',collapse=''),header=T,sep='\t',comment.char="",stringsAsFactors = F)
  }
  
  if (is.numeric(gene_of_interest)) {
    out_cnt=rep(0,nrow(out_mat[[1]]))
    for (i in 1:nrow(files)) {
      out_cnt=out_cnt+  apply(out_mat[[i]][,-1],1,sum)  }
    top_ind=order(out_cnt,decreasing = T)
    gene_of_interest=out_mat[[1]]$gene[top_ind[1:gene_of_interest]]
  }
  
  temp=data.frame(matrix(NaN,ncol(out_mat[[1]])-1,nrow(files)+2))
  temp[,1]=names(out_mat[[1]])[-1]
  ind=match(gene_of_interest,out_mat[[1]]$gene)
  for (i in 1:length(ind)) {
    cur_gene=out_mat[[1]][ind[i],1]
    filename=paste(img_dir,'/', direction,'_',cur_gene,'.png',sep='')
    if (!file.exists(filename)) {
    for (j in 1:nrow(files)) {
      temp[,j+2]=t(j*out_mat[[j]][-1][ind[i],]) }
    names(temp)=c('Subj','Subtype',files[,2])
    temp=temp[,c(1,ncol(temp):2)]
    temp[,ncol(temp)]=5
    color_levels=match(names(temp),files[,2])
    color_levels=color_levels[!is.na(color_levels)]
    color_levels=as.character(c(0,color_levels,5))
    df=melt(temp)
    df$value=as.factor(df$value)
    cols=c('0'="#D9D9D9",files[,4],"red")
    ggplot(df,aes(Subj,variable,fill=value)) + 
      geom_tile(stat='identity',aes(width=0.8,height=0.8))+
      scale_fill_manual(values=c("#D9D9D9",files[,4],"red"),limits=color_levels)+
      labs(title=paste(cur_gene,'_',direction,sep=''))+
      coord_fixed() + theme(legend.position="none")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_text(color='black',size=floor(15/nrow(files))),
            plot.title=element_text(hjust=0.5,vjust=2,size=10),
            panel.background = element_blank())
    ggsave(file=filename,device='png',dpi=800) }}}


enrich_outlier=function(suffix,direction) {
  
  files=read.table(paste('/input/data/sample_input_', suffix,'.txt',sep='',collapse=''),header=F,sep='\t',comment.char="",stringsAsFactors = F)
  clinical=read.table('/input/data/CPTAC Ovary Clinical Data_20160915.csv',header=T,sep=',',comment.char="",stringsAsFactors = F)
  clinical=clinical[clinical$Tumor.Grade !='' | clinical$Tumor.Stage..Pathological..Ovary.FIGO.Staging.System!='' ,]
  
  out_mat=list()
  for (i in 1:nrow(files)) {
    out_mat[[i]]=read.table(paste('/output/res/result/', suffix,'/', files[i,2],'_outlier_',direction, '.txt',sep='',collapse=''),header=T,sep='\t',comment.char="",stringsAsFactors = F)
    data=out_mat[[i]]
    names(data)=gsub('X','',names(data))
    ind= match(names(data),clinical$Participant.ID)
    
    phenotype1=clinical$Tumor.Stage..Pathological..Ovary.FIGO.Staging.System[ind]
    phenotype1[grep('Unknown',phenotype1,ignore.case=T)]=NA
    
    phenotype2=clinical$Tumor.Grade[ind]
    phenotype2[grep('Unknown|GX',phenotype2,ignore.case=T)]=NA
    phenotype=list(tumor_stage=phenotype1,tumor_grade=phenotype2)
    subtype1=list(statge1vs234=c('IC'),stage12vs34=c('IC','IIB'),stage123vs4=c('IC','IIB','III','IIIA','IIIB','IIIC'))
    subtype2=list(grade12vs3=c('G1','G2'),grade1vs23=c('G1'))
    subtype=list(subtype1,subtype2)
    result=matrix(data=NA,nrow=nrow(data),ncol=length(subtype[[1]])+length(subtype[[2]]))
    rownames(result)=data[,1]
    colnames(result)=c(names(subtype[[1]]),names(subtype[[2]]))
    cnt=0
    for (j in 1:length(phenotype)) {
      cur_phenotype=phenotype[[j]]
      for (k in 1:length(subtype[[j]])) {
        cnt=cnt+1
        cur_subtype=subtype[[j]][[k]]
        cat('\n',names(subtype[[j]])[k],'\n')
        
        for (m in 1:nrow(data)) {
          if (m%%100==0)
            cat('.')
          if (m%%3000==0)
            cat('.')
          cur_outlier=data[m,!(cur_phenotype %in% cur_subtype) & (!is.na(cur_phenotype)) & cur_phenotype !='']
          if (length(cur_outlier)==0)
            cur_outlier=-1
          other_outlier=data[m,cur_phenotype %in% cur_subtype]
          if (length(other_outlier)==0)
            other_outlier=-1
          tab=matrix(c(sum(cur_outlier>0),sum(cur_outlier==0),sum(other_outlier>0),sum(other_outlier==0)),2,2)
          rownames(tab)=c('outlier','not outlier')
          colnames(tab)=c('other',names(subtype[[j]])[k])
          tt=fisher.test(tab,alternative='greater')
          result[m,cnt]=tt$p.value
        }}
      
    }
    
    ind=apply(result<0.05,1,sum,na.rm=T)>0
    if (sum(ind)) {
      result2=data.frame(Genes=rownames(result),result)
      result=result2[ind,]
      res_dir=paste('/output/res/enrich_result/',suffix,sep='',collaspse='')
      filename=paste(res_dir,'/', files[i,2],'_enrich_',direction, '.txt',sep='',collapse='')
      
      if (!dir.exists(res_dir))
        dir.create(res_dir,recursive=T)
      if (!exists(filename))
        write.table(result,filename,sep='\t',row.names=F,col.names=T) } }}

options(stringsAsFactors=FALSE)

datanames=c('jhu','pnnl')
directions=c('up','down','both')
top=5
for (i in 1:length(datanames))
  for (j in 1:length(directions)) {
    get_outlier(datanames[i],directions[j])
    plot_outlier(datanames[i],directions[j],top)
    enrich_outlier(datanames[i],directions[j])}


