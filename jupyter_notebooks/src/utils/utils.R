
generate.PDF <- function(filename=filename,outputpath=outputpath,plot= plot,
                         width = 8.5,height = 5) {    
  pdf(paste0(outputpath,"/",filename), 
      width=width, height=height,onefile=T)
  print(plot)
  dev.off()
}

#Note: we get the minus cor! since the minus score can respresent the oxphos level
getCor<-function(mergedData,columnLoc=ncol(mergedData),method=c('spearman','pearson')){
  #remove the 1st and last col
  outcome=as.data.frame(matrix(NA,nrow = (ncol(mergedData)-2),ncol = 4))
  colnames(outcome)=c('object','Cor','pvalue','padj')
  outcome[,'object']=colnames(mergedData)[2:(ncol(mergedData)-1)]
  #remove 1st and last cols
  for (i in 1:(ncol(mergedData)-2)) {
    outcome[i,'Cor']=-cor(mergedData[,i+1],mergedData[,columnLoc],method = method)
    outcome[i,'pvalue']=cor.test(mergedData[,i+1],mergedData[,columnLoc],method = method)$p.value
  }
  outcome[,'padj']=p.adjust(outcome[,'pvalue'])
  return(outcome)
}
#run wilcox 
#Cor_FC_Expr: object and Group

options(stringsAsFactors = F)
plotCor<-function(data=cordata,x='',y='',title=''){
  require(ggplot2)
  data[,x]=as.numeric(data[,x])
  data[,y]=as.numeric(data[,y])
  
  data=data[!is.na(data[,x]),]
  tmp=cor.test(data[,x],data[,y],method = 'spearm')
  cor=tmp$estimate
  pvalue=tmp$p.value
  q=ggplot(data,aes_string(x=x,y=y))+geom_point()+geom_smooth(method = 'lm')+
    annotate(geom="text", x=max(data[,x])*4/5, y=max(data[,y])*4/5, label=paste('Correlation is',round(cor,3)))+
    annotate(geom="text", x=max(data[,x])*4/5, y=max(data[,y])*3/5, label=paste('P value is',round(pvalue,4)))+
    labs(y=y,x=x)+ggtitle(title)
  plot(q)
  
}

runwilcox<-function(Drugname=EgfrDrug,Groupinfo=Cor_FC_Expr,IC50=prob$mat,plot=T,
                    plotcor=T,threshold=0.01){
  require(ggpubr)
  Drug_IC50=IC50[rownames(IC50)%in% Drugname,]
  Drug_IC50=as.matrix(t(Drug_IC50))
  Drug_IC50_group=merge(Drug_IC50,Groupinfo,by.x=0,by.y='object')
  Drug_IC50_group$Group=as.factor(Drug_IC50_group$Group)
  OutMat=data.frame(drug=Drugname,pvalue=NA,FoldChange=NA)
  for (drug in Drugname) {
    tmp=wilcox.test(as.numeric(Drug_IC50_group[,drug])~ Drug_IC50_group$Group,data = Drug_IC50_group,paired = F)
    OutMat[OutMat$drug==drug,'pvalue']=tmp$p.value
    plotData=Drug_IC50_group[,c(drug,'Group')]
    colnames(plotData)[1]='IC50'
    plotData$IC50=as.numeric(plotData$IC50)
    #get FC
    tmp1=aggregate(.~Group,plotData,median)
    #FC=high-low
    OutMat[OutMat$drug==drug,'FoldChange']=tmp1$IC50[1]-tmp1$IC50[2]
    #plot significant outcome
    if (plot==T & tmp$p.value<threshold) {
      p <- ggboxplot(plotData, x = "Group", y = "IC50",
                     color = "Group", palette = "jco",add = "jitter")
      p<-p + stat_compare_means()+ggtitle(paste('Boxplot of IC50 level of ',drug,sep = ''))
      plot(p)
      #plot cor
      if (plotcor==T) {
        cordata=Drug_IC50_group[,c(drug,'Cor')]
        plotCor(cordata,x=drug,y='Cor')
      }
    }
  }
  OutMat=OutMat[order(OutMat$pvalue),]
  return(OutMat)
}

#geneexpr must have gene_id column and no other non-digital columns;Genesignature:Gene.symbol+FC is enough
estimateOxphos.ch<-function(GeneExpr,Genesignature,oxphos.prob=0.5,columnLoc='FC_treatment',normalize=T){
  GeneExpr=as.data.frame(GeneExpr)
  #estimate the OXPHOS level of each patient
  loc=which(GeneExpr$gene_id %in% Genesignature$Gene.symbol)
  GS_Expr_pd1=GeneExpr[loc,]
  rownames(GS_Expr_pd1)=GS_Expr_pd1[,'gene_id']
  GS_Expr_pd1=GS_Expr_pd1[,colnames(GS_Expr_pd1)!=c('gene_id')]
  #note!!!
  if (normalize==T) {
    GS_Expr_pd1=log2(GS_Expr_pd1+1)-rowMeans(log2(GS_Expr_pd1+1))
  }
  
  GS_Expr_pd1=GS_Expr_pd1[!is.na(GS_Expr_pd1[,1]),]
  GS_Expr_pd1_FC=merge(GS_Expr_pd1,Genesignature,by.x=0,by.y='Gene.symbol')
  loc=which(colnames(GS_Expr_pd1_FC)%in%c("average_DMSO_24_48","average_treatment"))
  if (length(loc)>0) {
    GS_Expr_pd1_FC=GS_Expr_pd1_FC[,-loc]
  }
  
  Cor_FC_Expr_pd1=getCor(mergedData = GS_Expr_pd1_FC,columnLoc = columnLoc,method = 'spearman')
  
  #classify patient into two groups
  Cor_FC_Expr_pd1$Group='Low'
  meidan.cor=quantile(Cor_FC_Expr_pd1$Cor,probs = oxphos.prob,na.rm = T)
  Cor_FC_Expr_pd1[Cor_FC_Expr_pd1$Cor> meidan.cor,]$Group='High'
  
  return(Cor_FC_Expr_pd1)
  
}

getSur_Res<-function(clinical_Cor,ciborsort,selCell,prob=0.25,response=T,cohort='Miao',
                     color = c("red","pink","darkgreen","green")){
  #plot surv curve and calc the response rate
  #clinical_cor :1st col is sample id,Group column ,sur+status(OS,OS.Event), Response is neccessary 
  #                 usually merged from estimate_oxphos results and clinical data
  #cibersort:1st col is "Input.Sample": sample id 
  require(ggplot2)
  require(survival)
  require(ggfortify)
  resposeRate=list()
  data.stat=list()
  clinical_Cor_cell.list=list()
  plot.list=list()
  for (cell in selCell) {
    
    #merge survive info
    clinical_Cor_cell=merge(clinical_Cor,ciborsort[,c("Input.Sample",cell)],by.x='object',by.y= "Input.Sample")
    clinical_Cor_cell$Group=paste('OXPHOS',clinical_Cor_cell$Group,sep = '_')
    clinical_Cor_cell$Cell=paste(cell,'High',sep = '_')
    #choose the criteria of classifying cell infiltration 
    cell.prob=quantile(clinical_Cor_cell[,cell],probs = prob)
    clinical_Cor_cell[clinical_Cor_cell[,cell]<= cell.prob, 'Cell']=paste(cell,'Low',sep = '_')
    clinical_Cor_cell$Group_2=paste(clinical_Cor_cell$Group,clinical_Cor_cell$Cell,sep = '+')
    #note~~~OS.event
    km_fit <- survfit(Surv(OS,event = OS.Event) ~ Group_2, data = clinical_Cor_cell)
    km_diff=survdiff(Surv(OS,event = OS.Event) ~ Group_2,data = clinical_Cor_cell)
    
    #Kaplan-Meier Curve :
    p=plotSur(fit =km_fit,fit.diff = km_diff, xlab = "overall survival (Month)",
              title = cohort,color = color)
    plot(p)
    plot.list[[cell]]= p
    if (response==T) {
      #calc responder's rate
      clinical_Cor_cell=clinical_Cor_cell[!is.na(clinical_Cor_cell$Response),]
      totalR=sum(clinical_Cor_cell$Response)
      totalNR=nrow(clinical_Cor_cell)-totalR
      r=aggregate(Response~Group_2,data = clinical_Cor_cell,sum)
      total=table(clinical_Cor_cell$Group_2)
      nr=total-r$Response
      data.resp=data.frame(pvalue=rep(NA,4))
      rownames(data.resp)=r$Group_2
      for (gp in r$Group_2) {
        tmp=rbind(c(r[r$Group_2==gp,2],nr[gp]),c(totalR-r[r$Group_2==gp,2],totalNR-nr[gp]))
        data.resp[gp,'pvalue']=fisher.test(tmp)$p.value
        data.stat[[paste(cell,gp,sep = '+')]]=tmp
      }
      resposeRate[[cell]]=data.resp
      
    }
    clinical_Cor_cell.list[[cell]]=clinical_Cor_cell
  }
  output=list(clinical_Cor_cell=clinical_Cor_cell.list,
              resposeRate=resposeRate,
              plot = plot.list)
  return(output)
}


#plot the DAVID results
plotDavidResults<-function(davidResult,num=10,title="Top pathway enrichment"){
  pathway=head(davidResult,n=num)
  colnames(pathway) = gsub(pattern = ' ',replacement = '_',
                         colnames(pathway))
  colnames(pathway)[4] = 'GeneRatio'
  pathway <- pathway %>% mutate(Term = str_remove(Term,pattern = 'GO:.*~')) %>%
    arrange(GeneRatio) 
  pathway$Term = factor(pathway$Term,
                        levels = pathway$Term)
  #points' size: the enriched gene number; points' color:-log10(FDR)
  q= ggplot(pathway,aes_string('Fold_Enrichment','Term'))+
    geom_point(aes(size=GeneRatio,color=-log10(Benjamini)))+
    scale_colour_gradient(low="blue",high="red")+
    labs(x="Fold Enrichment",y="Pathway name",title=title)+
    theme_bw(base_size = 14)
  plot(q)
}

DEgeneAnalysis<-function(TPM,Genesignature,countData,normalize=F,prob=0.1){
  require(data.table)
  require(DESeq2)
  
  TPM.oxphos=estimateOxphos.ch(GeneExpr =TPM,Genesignature = Genesignature,normalize = normalize )
  
  #redefine the group:topest 10 and lowest 10 cells
  print("Define oxphos group")
  top.threshold=quantile(TPM.oxphos$Cor,(1-prob))
  low.threshold=quantile(TPM.oxphos$Cor,prob)
  top.ccle=TPM.oxphos[TPM.oxphos$Cor>top.threshold,1,drop=F]
  top.ccle$group='Top'
  low.ccle=TPM.oxphos[TPM.oxphos$Cor<low.threshold,1,drop=F]
  low.ccle$group='Low'
  coldata.ccle=rbind(top.ccle,low.ccle)
  rownames(coldata.ccle)=coldata.ccle$object
  
  #DE analysis
  #1.import countData
  print(paste("Identify the DEgene between top and bottom ",prob,sep = " "))
  count.ccle.comp<- countData[, rownames(coldata.ccle)]
  all(rownames(coldata.ccle)== colnames(count.ccle.comp))
  #deseq2
  dds.ccle <- DESeqDataSetFromMatrix(countData =round(as.matrix(count.ccle.comp)),
                                     colData = coldata.ccle,
                                     design = ~ group)
  keep <- rowSums(counts(dds.ccle)) >= 15
  dds.ccle <- dds.ccle[keep,]
  #exec in R not rmarkerdown
  # parallel=TRUE for multiple thread
  print('processing deseq2...')
  dds.ccle <- DESeq(dds.ccle,parallel =T) 
  res <- results(dds.ccle, contrast=c("group","Top","Low"))
  resultsNames(dds.ccle)
  degene=as.data.frame(res)
  loc=which(degene$padj<0.01 & abs(degene$log2FoldChange)>1)
  degene=degene[loc,]
  out=list(Group_info=coldata.ccle,Degene=degene)
  return(out)
}

DEgeneAnalysis_GSVA<-function(GSVA_result,countData){
  require(data.table)
  require(DESeq2)
  #redefine the group:topest 10 and lowest 10 cells

  top.threshold=quantile(GSVA_result$ERRGeneSet,0.9)
  low.threshold=quantile(GSVA_result$ERRGeneSet,0.1)
  top=data.frame(object=rownames(GSVA_result[GSVA_result$ERRGeneSet>top.threshold,]),
                 group='Top')
  low=data.frame(object=rownames(GSVA_result[GSVA_result$ERRGeneSet<low.threshold,]),
                 group='Low')
  
  coldata.ccle=rbind(top,low)
  rownames(coldata.ccle)=coldata.ccle$object
  
  #DE analysis
  #1.import countData
  count.ccle.comp<- countData[, rownames(coldata.ccle)]
  all(rownames(coldata.ccle)== colnames(count.ccle.comp))
  #deseq2
  dds.ccle <- DESeqDataSetFromMatrix(countData =round(as.matrix(count.ccle.comp)),
                                     colData = coldata.ccle,
                                     design = ~ group)
  keep <- rowSums(counts(dds.ccle)) >= 15
  dds.ccle <- dds.ccle[keep,]
  #exec in R not rmarkerdown
  # parallel=TRUE for multiple thread
  print('processing deseq2...')
  dds.ccle <- DESeq(dds.ccle,parallel =T) 
  res <- results(dds.ccle, contrast=c("group","Top","Low"))
  resultsNames(dds.ccle)
  degene=as.data.frame(res)
  loc=which(degene$padj<0.01 & abs(degene$log2FoldChange)>1)
  degene=degene[loc,]
  out=list(Group_info=coldata.ccle,Degene=degene)
  return(out)
}

DEgeneAnalysis_PCA<-function(pca_score,countData){
  require(data.table)
  require(DESeq2)
  #redefine the group:topest 10 and lowest 10 cells
  
  top.threshold=quantile(pca_score$score,0.9)
  low.threshold=quantile(pca_score$score,0.1)
  top=data.frame(object=pca_score[pca_score$score>top.threshold,1],
                 group='Top')
  low=data.frame(object=pca_score[pca_score$score<low.threshold,1],
                 group='Low')
  
  coldata.ccle=rbind(top,low)
  rownames(coldata.ccle)=coldata.ccle$object
  
  #DE analysis
  #1.import countData
  count.ccle.comp<- countData[, rownames(coldata.ccle)]
  all(rownames(coldata.ccle)== colnames(count.ccle.comp))
  #deseq2
  dds.ccle <- DESeqDataSetFromMatrix(countData =round(as.matrix(count.ccle.comp)),
                                     colData = coldata.ccle,
                                     design = ~ group)
  keep <- rowSums(counts(dds.ccle)) >= 15
  dds.ccle <- dds.ccle[keep,]
  #exec in R not rmarkerdown
  # parallel=TRUE for multiple thread
  print('processing deseq2...')
  dds.ccle <- DESeq(dds.ccle,parallel =T) 
  res <- results(dds.ccle, contrast=c("group","Top","Low"))
  resultsNames(dds.ccle)
  degene=as.data.frame(res)
  loc=which(degene$padj<0.01 & abs(degene$log2FoldChange)>1)
  degene=degene[loc,]
  out=list(Group_info=coldata.ccle,Degene=degene)
  return(out)
}



compareMHC<-function(GeneExpr=skcm.tcga,Cor=skcm.tcga.cor,title='TCGA_SKCM',plot=T,Geneset=geneset){
  #rownames:Geneexpr
  loc=which(rownames(GeneExpr)%in% Geneset)
  
  #full word grep('^PDCD1$', rownames(skcm.tcga),value = T)
  HLAgenes_expr=GeneExpr[loc,]
  Geneset= rownames(HLAgenes_expr)
  HLAgenes_expr_cor=merge(t(HLAgenes_expr),Cor[,c('object','Cor','Group')],by.x=0,by.y=1)
  apply(HLAgenes_expr_cor[,colnames(HLAgenes_expr_cor)%in%Geneset], 2, as.numeric)->HLAgenes_expr_cor[,colnames(HLAgenes_expr_cor)%in%Geneset]
  
  #must change the gene name, or cant recongzed by ggboxplot
  Geneset=gsub(Geneset,pattern = '[-]',replacement = '_')
  colnames(HLAgenes_expr_cor)=gsub(colnames(HLAgenes_expr_cor),pattern = '[-]',replacement = '_')
  
  pvalueMat=data.frame(pvalue=rep(NA,length(Geneset)))
  rownames(pvalueMat)=Geneset
  for (name in Geneset) {
    high.group=HLAgenes_expr_cor[HLAgenes_expr_cor$Group=='High',name]
    low.group=HLAgenes_expr_cor[HLAgenes_expr_cor$Group=='Low',name]
    pvalue=wilcox.test(high.group,low.group)$p.value
    pvalueMat[name,'pvalue']=pvalue
    pvalueMat[name,'FC']=median(high.group)/median(low.group)
    pvalueMat=pvalueMat[!is.na(pvalueMat$FC),]
    pvalueMat=pvalueMat[pvalueMat$FC != Inf,]
    if(is.na(pvalue<0.05)) {next}
    if(pvalue<0.05&plot==T){
      require(ggpubr)
      p=ggboxplot(HLAgenes_expr_cor[,c(name,"Group")], x = "Group", y =name,color = 'Group',
                  add = "jitter",alpha=0.6)+ stat_compare_means()+ggtitle(title)
      plot(p)
      }
    
  }
  return(pvalueMat)
}



compareMHC.GSEA<-function(GeneExpr=skcm.tcga,GSVA=GSVA_result,title='TCGA_SKCM',plot=T,Geneset=geneset){
  #rownames:Geneexpr
  loc=which(rownames(GeneExpr)%in% Geneset)
  #full word grep('^PDCD1$', rownames(skcm.tcga),value = T)
  HLAgenes_expr=GeneExpr[loc,]
  Geneset= rownames(HLAgenes_expr)
  HLAgenes_expr_cor=merge(t(HLAgenes_expr),GSVA,by.x=0,by.y=0)
  apply(HLAgenes_expr_cor[,colnames(HLAgenes_expr_cor)%in%Geneset], 2, as.numeric)->HLAgenes_expr_cor[,colnames(HLAgenes_expr_cor)%in%Geneset]
  
  #must change the gene name, or cant recongzed by ggboxplot
  Geneset=gsub(Geneset,pattern = '[-]',replacement = '_')
  colnames(HLAgenes_expr_cor)=gsub(colnames(HLAgenes_expr_cor),pattern = '[-]',replacement = '_')
  
  pvalueMat=data.frame(pvalue=rep(NA,length(Geneset)))
  rownames(pvalueMat)=Geneset
  for (name in Geneset) {
    high.group=HLAgenes_expr_cor[HLAgenes_expr_cor$Group=='High',name]
    low.group=HLAgenes_expr_cor[HLAgenes_expr_cor$Group=='Low',name]
    pvalue=wilcox.test(high.group,low.group)$p.value
    pvalueMat[name,'pvalue']=pvalue
    pvalueMat[name,'change']='Low'
    if (median(high.group)>median(low.group)){
      pvalueMat[name,'change']='High'
    }
    
    if(is.na(pvalue<0.05)) {next}
    if(pvalue<0.05&plot==T){
      require(ggpubr)
      p=ggboxplot(HLAgenes_expr_cor[,c(name,"Group")], x = "Group", y =name,color = 'Group',
                  add = "jitter")+ stat_compare_means()+ggtitle(title)
      plot(p)
    }
    
  }
  return(pvalueMat)
}


#note: countdata and genename must match
estOxphosPCA<-function(countdata=countdata,Genesignature=genename,threshold=0.5,normalization=T){
  require(IMvigor210CoreBiologies)
  if (normalization==T) {
    #filter and normalized countdata
    voomD <- filterNvoom(countdata,minCpm=0.25)
    m <- voomD$E
    m <- t(scale( t( m ),
                  center=TRUE, 
                  scale=TRUE)
           )
  }
  else {m = countdata}
  tmp=m[rownames(m)%in%Genesignature,,drop=F]
  #only on value for the whole sg
  sg_pca_score <- gsScore(tmp)
  
  sg_pca_score=data.frame(object=names(sg_pca_score),score=sg_pca_score)
  sg_pca_score$Group='High'
  value=quantile(sg_pca_score$score,probs = threshold)
  sg_pca_score[sg_pca_score$score<value,'Group']='Low'
  return(sg_pca_score)
}





Count2TPM <- function(countMat, idType = "Ensembl", organism="GRCh38")
{
  if(organism=="GRCh38")
    ensembl = read.delim("/liulab/zzeng/TiSig/static/Annotation/hg38/hg38_mart_export.txt", 
                         check.names = FALSE)
  else
    ensembl = read.delim("/liulab/zzeng/TiSig/static/Annotation/mm10/mm10_mart_export.txt", check.names = FALSE)
  ensembl$Length <- abs(ensembl$`Gene end (bp)` - ensembl$`Gene start (bp)`)
  if(toupper(idType) == "ENSEMBL"){
    len <- ensembl[match(rownames(countMat),ensembl$`Gene stable ID`), "Length"]
    rownames(countMat) = ensembl[match(rownames(countMat),ensembl$`Gene stable ID`), "Gene name"]
  }else if(toupper(idType) == "SYMBOL")
    len <- ensembl[match(rownames(countMat),ensembl$`Gene name`), "Length"]
  else
    stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")
  
  na_idx = which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
    countMat = countMat[!is.na(len),]
    len = len[!is.na(len)]
  }
  tmp <- as.matrix(countMat / len)
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))
  TPM = TPM[!duplicated(rownames(TPM)),]
  return(TPM)
}

#function to plot star heatmap
#merged.tcga.hla is tidyverse,pvalue,
#factor for x and y
plot_starheatmap <- function(data = merged.tcga.hla,x="cancer.order",y="gene",
                             fill="coefficients",pvalue='pvalue',
                             color.xtext = cus_color,
                             palette = "RdYlBu"){
  data$stars <- cut(data[,pvalue], 
                    breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), 
                    label=c("***", "**", "*", "")) 
  p=ggplot(data, 
           aes_string(x=x,y=y,fill=fill)) +
    geom_tile() + 
    #scale_fill_gradient2(high= color[1], mid= color[2], low=color[3])+
    #scale_fill_distiller(palette = palette)+
    scale_fill_viridis(option="magma")+
    labs(y=NULL, x=NULL, fill=fill)  +
    geom_text(aes(label=stars), color="black", size=3)+
    theme_bw() 
  if (!is.na(color.xtext)) {
    p=p+theme(axis.text.x=element_text(angle = -45, hjust = 0,color= color.xtext))
  }else{
    p=p+theme(axis.text.x=element_text(angle = -45, hjust = 0))
  }
  plot(p)
}

returnfit<- function(data= sur.data,subset_value=c("ALT","WT")){
  out=list()
  fit <- survfit(Surv(survive, event=vital_status) ~Group2, 
                 data=data,
                 subset=(BRAF== subset_value))
  fit.diff = survdiff(Surv(survive, event=vital_status) ~Group2, 
                      data=data,
                      subset=(BRAF==subset_value))
  out[["fit"]]=fit
  out[["fit.diff"]]= fit.diff
  return(out)
}


plotSur <- function(fit=fit1,title='Kaplan-Meier Curve of non-mutated BRAF SKCM',
                    xlab = "overall survival(Day)",fit.diff = fit1.diff,
                    color = c("red","green")
                    ){
  pvalue = round(pchisq(fit.diff$chisq, length(fit.diff$n)-1, lower.tail = FALSE),4)
  p=autoplot(fit,conf.int = FALSE,surv.size =1,
             xlab = xlab, ylab = "Probability of survival")+
    ggtitle(title)+
    theme_classic()+
    theme(legend.position = c(0.8, 0.8))+
    annotate(geom = 'text',x=0.8*max(fit$time),y=0.7,
             label=c(paste(' logrank pvalue ',pvalue,sep = '')))+
    scale_color_manual(values = color)
  plot(p)
}

