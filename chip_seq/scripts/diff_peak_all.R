# load packages
suppressMessages(library(DESeq2))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

# read in command line params
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="The path to data containing counts for the peaks", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="The path to output file", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# read in data
df <- read.delim(opt$input, sep='\t')
samps <- colnames(df)[4:length(colnames(df))]
groups <- unique(str_replace(str_sub(samps, 1, -4), '_sgRNA[12]', ''))

# establish baseline group
#if ('WT' %in% groups) {
  #wt_group <- 'WT'
#}
if ('ROSA' %in% groups) {
  wt_group <- 'ROSA'
}


# perform differential analysis
for (group in groups) {
  if (group != wt_group){
    # get expression object
    if (group=="COP1_KO"){
      tmp_samps <- c("ROSA_R1", "ROSA_R2",
                     paste0(group, '_sgRNA1_R1'), paste0(group, '_sgRNA1_R2')) # paste0(group, '_sgRNA2_R1'), paste0(group, '_sgRNA2_R2'))
      expr <- as.matrix(df[,tmp_samps])
      SampleAnn <- data.frame(condition=c(0, 0, 1, 1), row.names=tmp_samps)
    } else if (group=="WT"){
      tmp_samps <- c('ROSA_R1', 'ROSA_R2',
                     paste0(group, '_R3'), paste0(group, '_R4')) # paste0(group, '_sgRNA2_R1'), paste0(group, '_sgRNA2_R2'))
      expr <- as.matrix(df[,tmp_samps])
      SampleAnn <- data.frame(condition=c(0, 0, 1, 1), row.names=tmp_samps)
    } else {
      tmp_samps <- c('WT_R3', 'WT_R4',
                     paste0(group, '_R1'), paste0(group, '_R2'))
      expr <- as.matrix(df[,tmp_samps])
      SampleAnn <- data.frame(condition=c(0, 0, 1, 1), row.names=tmp_samps)
    }
    obj = ExpressionSet(assayData = expr, phenoData = AnnotatedDataFrame(SampleAnn))

    # create design matrix
    design = stats::model.matrix(~1+condition, pData(obj))
    rownames(design) = sampleNames(obj)

    # run DESeq2
    dds <- DESeqDataSetFromMatrix(exprs(obj), colData = pData(obj), design = design)
    dds <- DESeq(dds)
    res <- results(dds)
    #res <- lfcShrink(dds, coef = ncol(design), quiet = TRUE)
    res$padj[is.na(res$padj)] <- 1
    res <- as.data.frame(res[, c("log2FoldChange", "baseMean", "stat", "pvalue", "padj")])
    colnames(res) <- c("log2FC", "baseMean", "stat", "pvalue", "padj")

    # annotate differential peaks
    res[group] <- 0
    is_signif <- res$padj<=0.05 & abs(res$log2FC)>0 #1
    res[is_signif & res['log2FC']< -0.5, group] <- -1
    res[is_signif & res['log2FC']> 0.5, group] <- 1
    df[group] <- res[group]
    df[paste(group, 'log2FC')] <- res$log2FC
    df[paste(group, "pvalue")] <- res$pvalue
    df[paste(group, 'padj')] <- res$padj
  }
}

# save results
write.table(df, opt$output, sep='\t', quote=F, row.names=F)

# figure out significant up peaks
#is_ge_up <- df['MUT_GE']==1
#is_am_up <- df['MUT_AM']==1
#is_400_up <- df['MUT_400']==1
#is_ge_down <- df['MUT_GE']==-1
#is_am_down <- df['MUT_AM']==-1
#is_400_down <- df['MUT_400']==-1


#samps <- c('WT_R1', 'WT_R2', 'MUT_AM_R1', 'MUT_AM_R2', 'MUT_GE_R1', 'MUT_GE_R2')
#expr <- as.matrix(df[,samps])
#SampleAnn <- data.frame(condition=c(1, 1, 0, 0), row.names=samps)


#obj = ExpressionSet(assayData = expr, phenoData = AnnotatedDataFrame(SampleAnn))

#tmp = GeneAnn[rownames(obj), drop = FALSE]
#rownames(tmp) = rownames(obj)
#slot(obj, "featureData") = AnnotatedDataFrame(tmp)


#design = stats::model.matrix(~1+condition, pData(obj))
#rownames(design) = sampleNames(obj)

# run DESeq2
#dds = DESeqDataSetFromMatrix(exprs(obj), colData = pData(obj), design = design)

#dds <- DESeq(dds)
#res <- lfcShrink(dds, coef = ncol(design), quiet = TRUE)
#res$padj[is.na(res$padj)] = 1
#res = res[, c("log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
#colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
