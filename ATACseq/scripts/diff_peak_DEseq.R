if(!require(optparse)){
  install.packages("optparse")
  suppressPackageStartupMessages(library(optparse))
}
# if(!require(DESeq2)){
#   install.packages("DESeq2")
#   suppressPackageStartupMessages(library(DESeq2))
# }
if(!require(DESeq2)){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install('DESeq2')
}

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(DESeq2))
option_list = list(
  make_option(c("-c", "--concount"), type="character", default=NULL, 
              help="Prefix of read count file of control samples, replicates should be seperate by comma.", metavar="character"),
  make_option(c("-t", "--treatcount"), type="character", default=NULL, 
              help="Prefix of read count file of treatment samples, replicates should be seperate by comma.", metavar="character"),
  make_option(c("-p", "--pvalue"), type="character", default="0.05", 
              help="Pvalue of diff peaks [default= %default]", metavar="character"),
  make_option(c("-f", "--logfoldchange"), type="character", default="1", 
              help="Log 2 fold change of diff peaks [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="diff_peak_out", 
              help="output file name [default= %default]", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#parse_args(OptionParser(option_list=option_list), args = c("--con_count=asd","--pvalue=0.05", "--logfoldchange=1"))
if (is.null(opt$concount)|is.null(opt$treatcount)){
  print_help(opt_parser)
  stop("At least one count files of each group must be supplied (input file).n", call.=FALSE)
}
# read count from htseq-count

con_file=opt$concount
treat_file=opt$treatcount
pvalue=as.numeric(opt$pvalue)
logfc=as.numeric(opt$logfoldchange)
out=opt$out
# tmp code
#con_file='Cop1sg1_IFN_r1_unique.sorted,Cop1sg1_IFN_r2_unique.sorted,Cop1sg2_IFN_r1_unique.sorted,Cop1sg2_IFN_r2_unique.sorted'
#treat_file='Rosa26_IFN_r1_unique.sorted,Rosa26_IFN_r2_unique.sorted'
#out='cop1_ifn_vs_rosa26_ifn'
# end

con_file_ls=unlist(strsplit(con_file,split = ','))
treat_file_ls=unlist(strsplit(treat_file,split = ','))
### input read count
con_count=data.frame()
n=1
for (c in con_file_ls){
  cf=read.table(paste0(out,'/count/',c,'.count.txt'),header = F,sep = '\t')
  cf=cf[-((nrow(cf)-4):nrow(cf)),]
  colnames(cf)=c('id',paste('con',n,sep = '_'))
  if (nrow(con_count)==0){
    con_count=rbind(con_count,cf)
  }else{
    con_count=merge(con_count,cf,by = 'id')
  }
  n=n+1
}

treat_count=data.frame()
n=1
for (c in treat_file_ls){
  tf=read.table(paste0(out,'/count/',c,'.count.txt'),header = F,sep = '\t')
  tf=tf[-((nrow(tf)-4):nrow(tf)),]
  colnames(tf)=c('id',paste('treat',n,sep = '_'))
  if (nrow(treat_count)==0){
    treat_count=rbind(treat_count,tf)
  }else{
    treat_count=merge(treat_count,tf,by = 'id')
  }
  n=n+1
}

### identify differential peaks using deseq2
all_count=merge(con_count,treat_count,by = 'id')
rownames(all_count)=all_count$id
all_count=all_count[,-1]

countData <- all_count
colData <- data.frame(row.names=colnames(countData),
                      condition=c(rep('con',(ncol(con_count)-1)),rep('treat',(ncol(treat_count)-1))),
                      libType= c(rep("single-end",2)))
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)

DESeq2.res=as.data.frame(res)
DESeq2.res.up=DESeq2.res[which(res$log2FoldChange>=logfc&DESeq2.res$pvalue<=pvalue),]
DESeq2.res.down=DESeq2.res[which(DESeq2.res$log2FoldChange<=(-logfc)&DESeq2.res$pvalue<=pvalue),]

if (!dir.exists(paste0(out,'/diff_peak/'))){
  dir.create(paste0(out,'/diff_peak/'))
}
write.table(DESeq2.res,file=paste0(out,'/diff_peak/DESeq2_res.txt'),quote = F,sep = '\t')
write.table(DESeq2.res.up,file=paste0(out,'/diff_peak/DESeq2_res_up.txt'),quote = F,sep = '\t')
write.table(DESeq2.res.down,file=paste0(out,'/diff_peak/DESeq2_res_down.txt'),quote = F,sep = '\t')

# ===== output peaks bed file
dd=read.table(paste0(out,'/count/merge_peak.gtf'))
plx_res_up_peak=dd[which(dd$V10%in%rownames(DESeq2.res.up)),]
plx_res_down_peak=dd[which(dd$V10%in%rownames(DESeq2.res.down)),]
plx_com_peak=dd[which(!dd$V10%in%c(rownames(DESeq2.res.down),rownames(DESeq2.res.up))),]


plx_res_up_peak=plx_res_up_peak[,c(-2,-3,-6,-7,-8,-9)]
plx_res_down_peak=plx_res_down_peak[,c(-2,-3,-6,-7,-8,-9)]
plx_com_peak=plx_com_peak[,c(-2,-3,-6,-7,-8,-9)]

plx_res_up_peak_center=plx_res_up_peak
plx_res_up_peak_center$start=format(round(rowMedians(as.matrix(plx_res_up_peak_center[,c(2,3)]))),scientific = F)
plx_res_up_peak_center$end=format(round(rowMedians(as.matrix(plx_res_up_peak_center[,c(2,3)]))+1),scientific = F)
plx_res_up_peak_center=plx_res_up_peak_center[,c(-2,-3,-4)]

plx_res_down_peak_center=plx_res_down_peak
plx_res_down_peak_center$start=format(round(rowMedians(as.matrix(plx_res_down_peak_center[,c(2,3)]))),scientific = F)
plx_res_down_peak_center$end=format(round(rowMedians(as.matrix(plx_res_down_peak_center[,c(2,3)]))+1),scientific = F)
plx_res_down_peak_center=plx_res_down_peak_center[,c(-2,-3,-4)]

plx_com_peak_center=plx_com_peak
plx_com_peak_center$start=format(round(rowMedians(as.matrix(plx_com_peak_center[,c(2,3)]))),scientific = F)
plx_com_peak_center$end=format(round(rowMedians(as.matrix(plx_com_peak_center[,c(2,3)]))+1),scientific = F)
plx_com_peak_center=plx_com_peak_center[,c(-2,-3,-4)]


write.table(plx_res_up_peak,paste0(out,'/diff_peak/up_peak.bed'),quote = F,sep = '\t',row.names = F,col.names = F)
write.table(plx_res_down_peak,paste0(out,'/diff_peak/down_peak.bed'),quote = F,sep = '\t',row.names = F,col.names = F)
write.table(plx_com_peak,paste0(out,'/diff_peak/com_peak.bed'),quote = F,sep = '\t',row.names = F,col.names = F)
write.table(plx_res_up_peak_center,paste0(out,'/diff_peak/up_peak_center.bed'),quote = F,sep = '\t',row.names = F,col.names = F)
write.table(plx_res_down_peak_center,paste0(out,'/diff_peak/down_peak_center.bed'),quote = F,sep = '\t',row.names = F,col.names = F)
write.table(plx_com_peak_center,paste0(out,'/diff_peak/com_peak_center.bed'),quote = F,sep = '\t',row.names = F,col.names = F)



