library(Seurat)
library(Signac)
library(tidyverse)
# library(ArchR)
library(future.apply)
library(ggpubr)
library(reshape2)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(Gviz)

fig_dir <- "figs/"

library(EnsDb.Hsapiens.v86)
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)
genebodyandpromoter.coords <- genebodyandpromoter.coords %>% subset(seqnames %in% c(1:22,'Y','X'))

linc_genes <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == 'lincRNA')

umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)


NucSeq.atac <- readRDS(file='mergeSamples_snATAC_TF_Gene_Activity_PartIII_Done.rds')
fig_dir <- "figures/"

cur_gene <- 'MBP'
cur_celltype <- 'ODC'
pdf(paste0(fig_dir, 'CoveragePlot_',  cur_gene, '.pdf'), width=4, height=5)
Signac::CoveragePlot(
  NucSeq.atac,
  region=cur_gene,
  group.by='monocle_clusters_umap_Cell.Type',
  extend.upstream=1000,
  extend.downstream=1000,
  peaks=FALSE
)
dev.off()



cur_pos <- subset(genebodyandpromoter.coords, symbol==cur_gene)
cur_position <- paste0('chr', as.character(seqnames(cur_pos)), ':',start(cur_pos)-2000,'-',start(cur_pos)) # promoter only
chr <- as.character(seqnames(cur_pos))
gen <- 'hg38'
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
print(paste(cur_gene, 'chr:', as.character(seqnames(cur_pos)), ',', abs(start(cur_pos) - end(cur_pos)) + 2000))

pdf(paste0(fig_dir, '/', cur_celltype, '/ideogram_', cur_gene, '.pdf'), width=8, height=2)
plotTracks(list(itrack), from=start(cur_pos), to=end(cur_pos), showId=FALSE)
dev.off()


