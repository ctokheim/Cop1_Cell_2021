suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
library(stringr)
library(optparse)

option_list = list(
  make_option(c("-m", "--deseq2_mat"), type="character", default=NULL, 
              help="signature reference", metavar="character"),
  make_option(c("-k", "--hallmark"), type="character", default=NULL, 
              help="msigdb hallmark file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

gsea_go <- function(data, outpath, type) {
    myGeneList <- data$stat
    names(myGeneList) <- data$ENTREZID
    myGeneList <- sort(myGeneList, decreasing = T)
    goEnrich <- gseGO(geneList     = myGeneList,
                       OrgDb        = org.Mm.eg.db,
                       keyType      = "ENTREZID",
                       ont          = type,
                       nPerm        = 100000,
                       pvalueCutoff = 1.05,
                       verbose      = FALSE
    )
    tmp <- as.data.frame(goEnrich)
    write.table(tmp, outpath, sep="\t", row.names=F)
}

gsea_kegg <- function(data, outpath) {
    myGeneList <- data$stat
    names(myGeneList) <- data$ENTREZID
    myGeneList <- sort(myGeneList, decreasing = T)
    keggEnrich <- gseKEGG(geneList     = myGeneList,
                        organism = "mmu",
                        nPerm        = 100000,
                        pvalueCutoff = 1.05,
                        verbose      = FALSE
    )
    tmp <- as.data.frame(keggEnrich)
    write.table(tmp, outpath, sep="\t", row.names=F)
}

gsea_hallmark <- function(data, hallmark, outpath) {
    myGeneList <- data$stat
    names(myGeneList) <- data$ENTREZID_human
    myGeneList <- sort(myGeneList, decreasing = T)
    hallmarkEnrich <- GSEA(myGeneList, TERM2GENE=h, verbose=FALSE, 
                           nPerm=100000, pvalueCutoff=1.05)
    tmp <- as.data.frame(hallmarkEnrich)
    write.table(tmp, outpath, sep="\t", row.names=F)
}

# read in data
df <- read.delim(opt$deseq2_mat)
tmp_gene_list <- df[,'gene_name']
geneToId <- bitr(tmp_gene_list, 
                 fromType="SYMBOL", 
                 toType="ENTREZID", 
                 OrgDb="org.Mm.eg.db")
# add in entrez id
df <- inner_join(df, geneToId, by=c('gene_name' = 'SYMBOL'))
df <- df[!is.na(df$ENTREZID) & !is.na(df$stat),]

# read in hallmark
h <- read.gmt(opt$hallmark)

# perform GSEA
bp_path <- paste0(opt$outdir, '/GO_BP.txt')
gsea_go(df, bp_path, type="BP")
mf_path <- paste0(opt$outdir, '/GO_MF.txt')
gsea_go(df, mf_path, type="MF")
cc_path <- paste0(opt$outdir, '/GO_CC.txt')
gsea_go(df, cc_path, type="CC")
kegg_path <- paste0(opt$outdir, '/KEGG.txt')
gsea_kegg(df, kegg_path)

# add human ID version
df['gene_name_human'] <- str_to_upper(df[,'gene_name'])
tmp_gene_list <- df[, 'gene_name_human']
geneToId <- bitr(tmp_gene_list, 
                 fromType="SYMBOL", 
                 toType="ENTREZID", 
                 OrgDb="org.Hs.eg.db")
df <- inner_join(df, geneToId, by=c('gene_name_human' = 'SYMBOL'), suffix=c('', '_human'))

# perform hallmark GSEA
hallmark_path <- paste0(opt$outdir, '/HALLMARK.txt')
gsea_hallmark(df, h, hallmark_path)
