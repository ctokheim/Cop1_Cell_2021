library(MAESTRO)
library(Seurat)
library(hdf5r)
library(Matrix)
library(ggplot2)
library(reshape)
library(cowplot)
library(RColorBrewer)

### B10
{
  matrix_dir = "../data/b10/cell_ranger/filtered_feature_bc_matrix/"
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  b10_mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(b10_mat) = barcode.names$V1
  rownames(b10_mat) = feature.names$V2
  
  
  ### separate samples with feature barcode
  
  feature_group=b10_mat[c("Rosa26_TotalSeqC1", "Cop1sg1_TotalSeqC2","Cop1sg2_TotalSeqC_3","unknown_TotalSeqC4"), ]
  feature_group=as.data.frame(t(as.data.frame(feature_group)))
  
  feature_group$group=colnames(feature_group)[max.col(feature_group,ties.method="first")]
  
  b10_Rosa26_barcode=rownames(feature_group)[which(feature_group$group=='Rosa26_TotalSeqC1')]
  b10_Cop1sg1_barcode=rownames(feature_group)[which(feature_group$group=='Cop1sg1_TotalSeqC2')]
  b10_Cop1sg2_barcode=rownames(feature_group)[which(feature_group$group=='Cop1sg2_TotalSeqC_3')]
  b10_unknown_barcode=rownames(feature_group)[which(feature_group$group=='unknown_TotalSeqC4')]

  rownames(b10_mat)=toupper(rownames(b10_mat))
}

### B11
{
  matrix_dir = "../data/b11/cellranger/filtered_feature_bc_matrix/"
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  b11_mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(b11_mat) = barcode.names$V1
  rownames(b11_mat) = feature.names$V2
  
  
  ### seperate samples with feature barcode
  
  feature_group=b11_mat[c("Reg_Rosa26_TotalSeqC1", "Reg_Cop1sg1_TotalSeqC2","Reg_Cop1sg2_TotalSeqC_3","Reg_unknown_TotalSeqC4"), ]
  feature_group=as.data.frame(t(as.data.frame(feature_group)))
  
  feature_group$group=colnames(feature_group)[max.col(feature_group,ties.method="first")]
  
  b11_Rosa26_barcode=rownames(feature_group)[which(feature_group$group=='Reg_Rosa26_TotalSeqC1')]
  b11_Cop1sg1_barcode=rownames(feature_group)[which(feature_group$group=='Reg_Cop1sg1_TotalSeqC2')]
  b11_Cop1sg2_barcode=rownames(feature_group)[which(feature_group$group=='Reg_Cop1sg2_TotalSeqC_3')]
  b11_unknown_barcode=rownames(feature_group)[which(feature_group$group=='Reg_unknown_TotalSeqC4')]
  
  rownames(b11_mat)=toupper(rownames(b11_mat))
  
}

### B12
{
  matrix_dir = "../data/b12/cellranger/filtered_feature_bc_matrix/"
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  b12_mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(b12_mat) = barcode.names$V1
  rownames(b12_mat) = feature.names$V2
  
  
  ### seperate samples with feature barcode
  
  feature_group=b12_mat[c("Pd1_Rosa26_TotalSeqC1", "Pd1_Cop1sg1_TotalSeqC2","Pd1_Cop1sg2_TotalSeqC_3","Pd1_unknown_TotalSeqC4"), ]
  feature_group=as.data.frame(t(as.data.frame(feature_group)))
  
  feature_group$group=colnames(feature_group)[max.col(feature_group,ties.method="first")]
  
  b12_Rosa26_barcode=rownames(feature_group)[which(feature_group$group=='Pd1_Rosa26_TotalSeqC1')]
  b12_Cop1sg1_barcode=rownames(feature_group)[which(feature_group$group=='Pd1_Cop1sg1_TotalSeqC2')]
  b12_Cop1sg2_barcode=rownames(feature_group)[which(feature_group$group=='Pd1_Cop1sg2_TotalSeqC_3')]
  b12_unknown_barcode=rownames(feature_group)[which(feature_group$group=='Pd1_unknown_TotalSeqC4')]
  
  rownames(b12_mat)=toupper(rownames(b12_mat))
}

### MAESTRO annotation cluster
data(human.immune.CIBERSORT)

### integrate cop1 and rosa26 with Seurat pipeline
{
  # mark rosa26 and cop1 ko
  # Initialize the Seurat object with the raw (non-normalized data).
  b10_mat_df=as.data.frame(b10_mat)
  colnames(b10_mat_df)=paste0('b10',colnames(b10_mat_df))
  b11_mat_df=as.data.frame(b11_mat)
  colnames(b11_mat_df)=paste0('b11',colnames(b11_mat_df))
  all_df<-merge(b10_mat_df,b11_mat_df,by=0)
  rownames(all_df)=all_df$Row.names;all_df=all_df[,-1]
  rm(b10_mat_df); rm(b11_mat_df)
  b12_mat_df=as.data.frame(b12_mat)
  colnames(b12_mat_df)=paste0('b12',colnames(b12_mat_df))
  all_df<-merge(all_df,b12_mat_df,by=0)
  rownames(all_df)=all_df$Row.names;all_df=all_df[,-1]
  rm(b12_mat_df)
  
  all_mat <- as(Matrix(as.matrix(all_df)), "dgTMatrix")
  
  rm(all_df)
  
  SeuratObj <- CreateSeuratObject(counts = all_mat, project = "all_intergrate", min.cells = 3, min.features = 200)
  SeuratObj
  ident=data.frame(row.names=c(paste0('b10',b10_Rosa26_barcode),
                               paste0('b10',b10_Cop1sg1_barcode),
                               paste0('b10',b10_Cop1sg2_barcode),
                               paste0('b11',b11_Rosa26_barcode),
                               paste0('b11',b11_Cop1sg1_barcode),
                               paste0('b11',b11_Cop1sg2_barcode),
                               paste0('b12',b12_Rosa26_barcode),
                               paste0('b12',b12_Cop1sg1_barcode),
                               paste0('b12',b12_Cop1sg2_barcode)),
                   ident = c(rep('B10_Rosa26',length(b10_Rosa26_barcode)),
                             rep('B10_Cop1sg1',length(b10_Cop1sg1_barcode)),
                             rep('B10_Cop1sg2',length(b10_Cop1sg2_barcode)),
                             rep('B11_Rosa26',length(b11_Rosa26_barcode)),
                             rep('B11_Cop1sg1',length(b11_Cop1sg1_barcode)),
                             rep('B11_Cop1sg2',length(b11_Cop1sg2_barcode)),
                             rep('B12_Rosa26',length(b12_Rosa26_barcode)),
                             rep('B12_Cop1sg1',length(b12_Cop1sg1_barcode)),
                             rep('B12_Cop1sg2',length(b12_Cop1sg2_barcode))))
  
  tmp3=SeuratObj@meta.data
  tmp3=merge(ident,tmp3,by = 0)
  rownames(tmp3)=tmp3$Row.names
  tmp3=tmp3[,c(-1,-3)]
  colnames(tmp3)[1]='orig.ident'
  SeuratObj@meta.data=tmp3
  
  # integration analysis
  SeuratObj.list <- SplitObject(SeuratObj, split.by = "orig.ident")
  
  SeuratObj.list <- lapply(X = SeuratObj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  immune.anchors <- FindIntegrationAnchors(object.list = SeuratObj.list, dims = 1:20)
  immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
  
  DefaultAssay(immune.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
  # t-SNE and Clustering
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
  immune.combined <- FindClusters(immune.combined, resolution = 0.5)
  
  immune.combined@project.name <- 'all_intergrate'
  
  ### Annotate each cluster
  # !!! time consuming step 
  set.seed(7)
  cluster.genes <- FindAllMarkersMAESTRO(object = immune.combined, 
                                         min.pct = 0.1, test.use = 'wilcox')
  saveRDS(cluster.genes,file = '../results/cluster.genes.rds')
  cluster.genes=readRDS('../results/cluster.genes.rds')
  # !!! time consuming step end
  
  cluster.genes <- cluster.genes[cluster.genes$p_val_adj < 1e-05, ]
  
  
  ann.immune.combined=list(RNA=immune.combined,gene=cluster.genes)
  
  ann.immune.combined$RNA <- RNAAnnotateCelltype(RNA = immune.combined, 
                                                 gene = ann.immune.combined$gene,
                                                 signatures = human.immune.CIBERSORT, 
                                                 min.score = 0.1)
  # UMAP plot
  p = DimPlot(object = ann.immune.combined$RNA, label = TRUE, pt.size = 0.2, 
              group.by = "assign.ident", label.size = 3, repel = T,cols = c( "#FB9A99" , "#FDBF6F", "#1F78B4" ,"#E31A1C" ,"#33A02C", "#FF7F00", "#CAB2D6" ,"#6A3D9A", "#B15928" , "#FF00FF", "#006666"))
  p
  #c("#A6CEE3" ,"#1F78B4" ,"#B2DF8A", "#33A02C", "#FB9A99" ,"#E31A1C" ,"#FDBF6F", "#FF7F00", "#CAB2D6" ,"#6A3D9A", "#FFFF99", "#B15928" ,"#E00F53", "#FF00FF", "#006666")
  ggsave(file.path(paste0('../results/',ann.immune.combined$RNA@project.name, "_group_ann.png")),p,width = 6.5,height = 4)

  ### merge cop1sg1 and cop1sg2 
  ann.immune.combined$RNA@meta.data$orig.ident_v2=gsub('sg[1-2]','',ann.immune.combined$RNA@meta.data$orig.ident)
  ### save combine object
  
  saveRDS(ann.immune.combined,file = '../results/ann.immune.combined.rds')
  ann.immune.combined=readRDS('../results/ann.immune.combined.rds')
  
}

### Visualize the change of the ratio of each cell type before and after Cop1 KO
{
  # merge sg1 and sg2
  anno_data=as.data.frame(ann.immune.combined$RNA$assign.ident)
  grou_data=as.data.frame(ann.immune.combined$RNA$orig.ident_v2)
  anno_data=merge(anno_data,grou_data,by=0)
  colnames(anno_data)=c('barcode','type','group')
  
  sum_fre_merge=data.frame()
  for (i in unique(anno_data$type)){
    for (j in unique(anno_data$group)){
      cell_num=nrow(anno_data[which(anno_data$group==j),])
      tmp=anno_data[which(anno_data$type==i&anno_data$group==j),]
      tmp_fre=data.frame(cell_type=i,cell_number=nrow(tmp),total_cells=cell_num,fre=nrow(tmp)/cell_num,group=j,channel=substr(j,1,3))
      sum_fre_merge=rbind(sum_fre_merge,tmp_fre)
    }
  }
  
  sum_fre_merge$group2=gsub('B*.+\\_','',sum_fre_merge$group)
  sum_fre_merge$group2=factor(sum_fre_merge$group2,levels = unique(sum_fre_merge$group2))
  dd_bar=sum_fre_merge
  #dd_bar$variable=factor(dd_bar$variable,levels = c('Rosa26','Cop1sg1','Cop1sg2','Cop1all'))
  # Stacked barplot with multiple groups
  p=ggplot(data=dd_bar, aes(x=channel, y=fre, fill=group2)) +
    geom_bar(stat="identity",position=position_dodge(),color='black')+theme_classic()+xlab('')+ylab('Frequency')+
    theme(legend.title=element_blank())+facet_wrap(~cell_type,ncol=2,scales = "free")+
    scale_fill_manual(values=c("#1F78B4" ,"#E31A1C" ,"#33A02C"))
    #scale_fill_brewer(palette="Dark2")
  p
  ggsave(filename = '../results/merge_s1_s2_fre_bar_merge.pdf',p,width = 8,height = 12)
}

