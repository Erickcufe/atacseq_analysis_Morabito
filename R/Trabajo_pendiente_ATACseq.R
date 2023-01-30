NucSeq.atac <- readRDS("ATAC_seurat_obj/mergeSamples_snATAC_TF_Gene_Activity_PartIII_Done.rds")
NucSeq.rna <- readRDS("ATAC_seurat_obj/datos_integrados_Anotados.rds")
NucSeq.rna$Cell.Type <- Idents(NucSeq.rna)
NucSeq.rna$tech <- 'rna'; NucSeq.atac$tech <- 'atac';
DefaultAssay(NucSeq.atac) <- 'RNA'

ft <- FindVariableFeatures(NucSeq.rna)
# saveRDS(NucSeq.atac, "mergeSamples_snATAC_TF_Gene_Activity_PartIII_Done.rds")
# compute anchors between RNA and ATAC
transfer.anchors <- FindTransferAnchors(
  reference=NucSeq.rna,
  query=NucSeq.atac,
  features=VariableFeatures(ft),
  reference.assay="RNA",
  query.assay="RNA",
  reduction="cca",
  verbose=T,
  dims=1:40
)
celltype.predictions <- TransferData(
  anchorset=transfer.anchors,
  refdata=NucSeq.rna$Cell.Type,
  weight.reduction=NucSeq.atac[["lsi"]],
  dims=1:40
)

NucSeq.atac <- AddMetaData(NucSeq.atac, celltype.predictions)
saveRDS(NucSeq.atac, "ATAC_seurat_obj/mergeSamples_cellTransfer.rds")

NucSeq.atac$predicted.id <- factor(NucSeq.atac$predicted.id, levels = levels(as.factor(NucSeq.rna$Cell.Type)))
Idents(NucSeq.atac) <- NucSeq.atac$predicted.id
Idents(NucSeq.rna) <- NucSeq.rna$Cell.Type
Idents(NucSeq.atac) <- paste0(Idents(NucSeq.atac), 'atac')
Idents(NucSeq.rna) <- paste0(Idents(NucSeq.rna), 'rna')
genes.use <- VariableFeatures(NucSeq.rna)
refdata <- GetAssayData(NucSeq.rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = NucSeq.atac[["lsi"]], dims = 1:50)

# Aqui nos quedamos 
NucSeq.atac[["RNA"]] <- imputation
NucSeq.atac <- RenameCells(NucSeq.atac, add.cell.id='atac')
NucSeq.rna <- RenameCells(NucSeq.rna, add.cell.id='rna')
NucSeq.coembed <- merge(x = NucSeq.rna, y = NucSeq.atac)
cells_to_keep <- c(colnames(NucSeq.rna), colnames(NucSeq.atac)[NucSeq.atac$prediction.score.max >= 0.5])
NucSeq.coembed <- NucSeq.coembed[,cells_to_keep]
NucSeq.coembed <- ScaleData(NucSeq.coembed, features = genes.use, do.scale = FALSE)
NucSeq.coembed <- RunPCA(NucSeq.coembed, features = rownames(NucSeq.coembed), verbose = FALSE)
NucSeq.coembed <- RunUMAP(NucSeq.coembed, dims = 1:30)
expr_matrix <- GetAssayData(NucSeq.coembed, slot='data', assay='RNA')
genes <- data.frame(as.character(rownames(expr_matrix)))
rownames(genes) <- rownames(expr_matrix)
genes <- as.data.frame(cbind(genes,genes))
colnames(genes) <- c("GeneSymbol", "gene_short_name")
NucSeq_cds <- new_cell_data_set(
  expr_matrix,
  cell_metadata=NucSeq.coembed@meta.data,
  gene_metadata=genes
)
NucSeq_cds@reducedDims[['PCA']] <- NucSeq.coembed@reductions$pca@cell.embeddings
NucSeq_cds <- align_cds(NucSeq_cds, preprocess_method='PCA', alignment_group = "Batch")
NucSeq_cds <- reduce_dimension(NucSeq_cds, reduction_method = 'UMAP', preprocess_method = "Aligned")
NucSeq_cds <- cluster_cells(NucSeq_cds, reduction_method='UMAP')



