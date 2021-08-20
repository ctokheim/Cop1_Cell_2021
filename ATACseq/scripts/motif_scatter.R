# if(!require(easyGgplot2)|!require(kassambara)){
#   library(devtools)
#   install_github("easyGgplot2", "kassambara")
#   suppressPackageStartupMessages(library(easyGgplot2))
# }else{
#   suppressPackageStartupMessages(library(easyGgplot2))
#   suppressPackageStartupMessages(library(ggrepel))
#   ggrepel
# }
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))
#library(devtools)
#install_github("easyGgplot2", "kassambara")

option_list = list(
  make_option(c("-f", "--motif_res"), type="character", default=NULL, 
              help="Motif enntichment result of Homer", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="Prefix of plot", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="motif", 
              help="output file name [default= %default]", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#parse_args(OptionParser(option_list=option_list), args = c("--con_count=asd","--pvalue=0.05", "--logfoldchange=1"))
if (is.null(opt$motif_res)){
  print_help(opt_parser)
  stop("Motif enriched file must be supplied (input file).n", call.=FALSE)
}

homer_file=opt$motif_res
prefix=opt$prefix
out=opt$out
# creat output directory
if (!dir.exists(paste0(out,'/homer/plot/'))){
  dir.create(paste0(out,'/homer/plot/'))
}

# input homer res

homer_res=read.table(homer_file,sep = '\t',header = F,skip = 1)
colnames(homer_res)=c('Motif_Name','Consensus','Pvalue','Log_Pvalue','q_value_Benjamini','Target_Sequences_with_Motif','Ratio_Target_Sequences_with_Motif','Background_Sequences_with_Motif','Ratio_Background_Sequences_with_Motif')

dd3=homer_res
colnames(dd3)[c(7,9)]=c('opt','ept')
dd3$opt=as.numeric(gsub('\\%','',dd3$opt))
dd3$ept=as.numeric(gsub('\\%','',dd3$ept))
dd3$neg_logP=(-dd3$Log_Pvalue)
dd3$genes=toupper(gsub('\\(+.*','',dd3$Motif_Name))
### scatterplot
p <- ggplot(dd3, aes(x=ept, y=opt, size=neg_logP)) +
  geom_point() + 
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()+
  geom_abline(slope = 1,intercept = 1,color='red',linetype='dashed')
p

p=p+xlab('Expected percentage of peaks')+ylab('Observed percentage of peaks')+
  theme (axis.text.x=element_text(angle = 0,hjust=1,size = 14,face = "bold"),
           axis.text.y=element_text(size = 14,face = "bold"),
           title=element_text(size = 14,face = "bold"),
           axis.title.x=element_text(size = 14,face = "bold"),
           legend.text=element_text(size = 14,face = "bold"),
           legend.title=element_text(size = 14,face = "bold"))
p

### select significant enriched motifs
dd3=dd3[order(dd3$Pvalue,decreasing = F),]
dd3_sig=dd3[which(dd3$q_value_Benjamini<=0.05),]

if (nrow(dd3_sig)>10){
  dat = dd3_sig[1:10,]
}else {
  dat = dd3_sig
}

if (nrow(dat)!=0){
  dat$group='sig'
  p <- p+geom_label_repel(
    data = dat,
    aes(label = genes,fill = group), color = 'black',
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
}
p

ggsave(plot = p,filename = paste0(out,'/homer/plot/',prefix,'.pdf'),width = 8,height = 6)


