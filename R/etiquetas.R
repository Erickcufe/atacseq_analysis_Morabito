library(Seurat)

so <- readRDS("SingleR_anotated_ATAC.rds")
Idents(so) <- so$Proplabels
DefaultAssay(so) <- 'ATAC'
DimPlot(so)

FeaturePlot(
  object = so,
  features = c('VIP', 'RORB', 'CDH9', 'SST', 'TRB1', 'SMAE3E'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

so_neurons <- so[,so$Proplabels == "Ex"]
