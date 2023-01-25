library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multicore", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

# Read peak of each sample of Morabito et al., 2021
peak_100 <- read.table(
  file = "../../samples_atacseq/sample-100-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_101 <- read.table(
  file = "../../samples_atacseq/sample-101-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_17 <- read.table(
  file = "../../samples_atacseq/sample-17-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_19 <- read.table(
  file = "../../samples_atacseq/sample-19-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_22 <- read.table(
  file = "../../samples_atacseq/sample-22-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_27 <- read.table(
  file = "../../samples_atacseq/sample-27-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_33 <- read.table(
  file = "../../samples_atacseq/sample-33-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_37 <- read.table(
  file = "../../samples_atacseq/sample-37-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_40 <- read.table(
  file = "../../samples_atacseq/sample-40-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_43 <- read.table(
  file = "../../samples_atacseq/sample-43-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_45 <- read.table(
  file = "../../samples_atacseq/sample-45-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_46 <- read.table(
  file = "../../samples_atacseq/sample-46-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_47 <- read.table(
  file = "../../samples_atacseq/sample-47-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_50 <- read.table(
  file = "../../samples_atacseq/sample-50-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_52 <- read.table(
  file = "../../samples_atacseq/sample-52-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_58 <- read.table(
  file = "../../samples_atacseq/sample-58-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_66 <- read.table(
  file = "../../samples_atacseq/sample-66-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_82 <- read.table(
  file = "../../samples_atacseq/sample-82-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_90 <- read.table(
  file = "../../samples_atacseq/sample-90-peaks.bed",
  col.names = c("chr", "start", "end")
)
peak_96 <- read.table(
  file = "../../samples_atacseq/sample-96-peaks.bed",
  col.names = c("chr", "start", "end")
)

# Convert to genomic ranges
gr.100 <- makeGRangesFromDataFrame(peak_100)
gr.101 <- makeGRangesFromDataFrame(peak_101)
gr.17 <- makeGRangesFromDataFrame(peak_17)
gr.19 <- makeGRangesFromDataFrame(peak_19)
gr.22 <- makeGRangesFromDataFrame(peak_22)
gr.27 <- makeGRangesFromDataFrame(peak_27)
gr.33 <- makeGRangesFromDataFrame(peak_33)
gr.37 <- makeGRangesFromDataFrame(peak_37)
gr.40 <- makeGRangesFromDataFrame(peak_40)
gr.43 <- makeGRangesFromDataFrame(peak_43)
gr.45 <- makeGRangesFromDataFrame(peak_45)
gr.46 <- makeGRangesFromDataFrame(peak_46)
gr.47 <- makeGRangesFromDataFrame(peak_47)
gr.50 <- makeGRangesFromDataFrame(peak_50)
gr.52 <- makeGRangesFromDataFrame(peak_52)
gr.58 <- makeGRangesFromDataFrame(peak_58)
gr.66 <- makeGRangesFromDataFrame(peak_66)
gr.82 <- makeGRangesFromDataFrame(peak_82)
gr.90 <- makeGRangesFromDataFrame(peak_90)
gr.96 <- makeGRangesFromDataFrame(peak_96)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.100, gr.101, gr.17, gr.19, gr.22, gr.27, gr.33, gr.37,
                               gr.40, gr.43, gr.45, gr.46, gr.47, gr.50, gr.52, gr.52,
                               gr.58, gr.66, gr.82, gr.90, gr.96))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# Load metadata
md.100 <- read.table(
  file = "barcodes_atacseq_AD/sample-100-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.101 <- read.table(
  file = "barcodes_atacseq_AD/sample-101-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.17 <- read.table(
  file = "barcodes_atacseq_AD/sample-17-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.19 <- read.table(
  file = "barcodes_atacseq_AD/sample-19-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.22 <- read.table(
  file = "barcodes_atacseq_AD/sample-22-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.27 <- read.table(
  file = "barcodes_atacseq_AD/sample-27-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.33 <- read.table(
  file = "barcodes_atacseq_AD/sample-33-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.37 <- read.table(
  file = "barcodes_atacseq_AD/sample-37-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.40 <- read.table(
  file = "barcodes_atacseq_AD/sample-40-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.43 <- read.table(
  file = "barcodes_atacseq_AD/sample-43-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.45 <- read.table(
  file = "barcodes_atacseq_AD/sample-45-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.46 <- read.table(
  file = "barcodes_atacseq_AD/sample-46-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.47 <- read.table(
  file = "barcodes_atacseq_AD/sample-47-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.50 <- read.table(
  file = "barcodes_atacseq_AD/sample-50-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.52 <- read.table(
  file = "barcodes_atacseq_AD/sample-52-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.58 <- read.table(
  file = "barcodes_atacseq_AD/sample-58-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.66 <- read.table(
  file = "barcodes_atacseq_AD/sample-66-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.82 <- read.table(
  file = "barcodes_atacseq_AD/sample-82-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.90 <- read.table(
  file = "barcodes_atacseq_AD/sample-90-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)
md.96 <- read.table(
  file = "barcodes_atacseq_AD/sample-96-singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)

# Perform an initial filtering of low count cells
md.100 <- md.100[md.100$passed_filters > 500, ]
md.101 <- md.101[md.101$passed_filters > 500, ]
md.17 <- md.17[md.17$passed_filters > 500, ]
md.19 <- md.19[md.19$passed_filters > 500, ]
md.22 <- md.22[md.22$passed_filters > 500, ]
md.27 <- md.27[md.27$passed_filters > 500, ]
md.33 <- md.33[md.33$passed_filters > 500, ]
md.37 <- md.37[md.37$passed_filters > 500, ]
md.40 <- md.40[md.40$passed_filters > 500, ]
md.43 <- md.43[md.43$passed_filters > 500, ]
md.45 <- md.45[md.45$passed_filters > 500, ]
md.46 <- md.46[md.46$passed_filters > 500, ]
md.47 <- md.47[md.47$passed_filters > 500, ]
md.50 <- md.50[md.50$passed_filters > 500, ]
md.52 <- md.52[md.52$passed_filters > 500, ]
md.58 <- md.58[md.58$passed_filters > 500, ]
md.66 <- md.66[md.66$passed_filters > 500, ]
md.82 <- md.82[md.82$passed_filters > 500, ]
md.90 <- md.90[md.90$passed_filters > 500, ]
md.96 <- md.96[md.96$passed_filters > 500, ]

# Create a fragment object
frags.100 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-100-fragments.tsv.gz",
  cells = rownames(md.100)
)
frags.101 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-101-fragments.tsv.gz",
  cells = rownames(md.101)
)
frags.17 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-17-fragments.tsv.gz",
  cells = rownames(md.17)
)
frags.19 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-19-fragments.tsv.gz",
  cells = rownames(md.19)
)
frags.22 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-22-fragments.tsv.gz",
  cells = rownames(md.22)
)
frags.27 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-27-fragments.tsv.gz",
  cells = rownames(md.27)
)
frags.33 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-33-fragments.tsv.gz",
  cells = rownames(md.33)
)
frags.37 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-37-fragments.tsv.gz",
  cells = rownames(md.37)
)
frags.40 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-40-fragments.tsv.gz",
  cells = rownames(md.40)
)
frags.43 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-43-fragments.tsv.gz",
  cells = rownames(md.43)
)
frags.45 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-45-fragments.tsv.gz",
  cells = rownames(md.45)
)
frags.46 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-46-fragments.tsv.gz",
  cells = rownames(md.46)
)
frags.47 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-47-fragments.tsv.gz",
  cells = rownames(md.47)
)
frags.50 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-50-fragments.tsv.gz",
  cells = rownames(md.50)
)
frags.52 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-52-fragments.tsv.gz",
  cells = rownames(md.52)
)
frags.58 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-58-fragments.tsv.gz",
  cells = rownames(md.58)
)
frags.66 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-66-fragments.tsv.gz",
  cells = rownames(md.66)
)
frags.82 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-82-fragments.tsv.gz",
  cells = rownames(md.82)
)
frags.90 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-90-fragments.tsv.gz",
  cells = rownames(md.90)
)
frags.96 <- CreateFragmentObject(
  path = "../../samples_atacseq/sample-96-fragments.tsv.gz",
  cells = rownames(md.96)
)

# Quantify peaks in each dataset
s100.counts <- FeatureMatrix(
  fragments = frags.100,
  features = combined.peaks,
  cells = rownames(md.100)
)
s101.counts <- FeatureMatrix(
  fragments = frags.101,
  features = combined.peaks,
  cells = rownames(md.101)
)
s17.counts <- FeatureMatrix(
  fragments = frags.17,
  features = combined.peaks,
  cells = rownames(md.17)
)
s19.counts <- FeatureMatrix(
  fragments = frags.19,
  features = combined.peaks,
  cells = rownames(md.19)
)
s22.counts <- FeatureMatrix(
  fragments = frags.22,
  features = combined.peaks,
  cells = rownames(md.22)
)
s27.counts <- FeatureMatrix(
  fragments = frags.27,
  features = combined.peaks,
  cells = rownames(md.27)
)
s33.counts <- FeatureMatrix(
  fragments = frags.33,
  features = combined.peaks,
  cells = rownames(md.33)
)
s37.counts <- FeatureMatrix(
  fragments = frags.37,
  features = combined.peaks,
  cells = rownames(md.37)
)
s40.counts <- FeatureMatrix(
  fragments = frags.40,
  features = combined.peaks,
  cells = rownames(md.40)
)
s43.counts <- FeatureMatrix(
  fragments = frags.43,
  features = combined.peaks,
  cells = rownames(md.43)
)
s45.counts <- FeatureMatrix(
  fragments = frags.45,
  features = combined.peaks,
  cells = rownames(md.45)
)
s46.counts <- FeatureMatrix(
  fragments = frags.46,
  features = combined.peaks,
  cells = rownames(md.46)
)
s47.counts <- FeatureMatrix(
  fragments = frags.47,
  features = combined.peaks,
  cells = rownames(md.47)
)
s50.counts <- FeatureMatrix(
  fragments = frags.50,
  features = combined.peaks,
  cells = rownames(md.50)
)
s52.counts <- FeatureMatrix(
  fragments = frags.52,
  features = combined.peaks,
  cells = rownames(md.52)
)
s58.counts <- FeatureMatrix(
  fragments = frags.58,
  features = combined.peaks,
  cells = rownames(md.58)
)
s66.counts <- FeatureMatrix(
  fragments = frags.66,
  features = combined.peaks,
  cells = rownames(md.66)
)
s82.counts <- FeatureMatrix(
  fragments = frags.82,
  features = combined.peaks,
  cells = rownames(md.82)
)
s90.counts <- FeatureMatrix(
  fragments = frags.90,
  features = combined.peaks,
  cells = rownames(md.90)
)
s96.counts <- FeatureMatrix(
  fragments = frags.96,
  features = combined.peaks,
  cells = rownames(md.96)
)

# Create Seurat Objects
s100_assay <- CreateChromatinAssay(s100.counts, fragments = frags.100)
s100 <- CreateSeuratObject(s100_assay, assay = "ATAC", meta.data = md.100)

s101_assay <- CreateChromatinAssay(s101.counts, fragments = frags.101)
s101 <- CreateSeuratObject(s101_assay, assay = "ATAC", meta.data = md.101)

s17_assay <- CreateChromatinAssay(s17.counts, fragments = frags.17)
s17 <- CreateSeuratObject(s17_assay, assay = "ATAC", meta.data = md.17)

s19_assay <- CreateChromatinAssay(s19.counts, fragments = frags.19)
s19 <- CreateSeuratObject(s19_assay, assay = "ATAC", meta.data = md.19)

s22_assay <- CreateChromatinAssay(s22.counts, fragments = frags.22)
s22 <- CreateSeuratObject(s22_assay, assay = "ATAC", meta.data = md.22)

s27_assay <- CreateChromatinAssay(s27.counts, fragments = frags.27)
s27 <- CreateSeuratObject(s27_assay, assay = "ATAC", meta.data = md.27)

s33_assay <- CreateChromatinAssay(s33.counts, fragments = frags.33)
s33 <- CreateSeuratObject(s33_assay, assay = "ATAC", meta.data = md.33)

s37_assay <- CreateChromatinAssay(s37.counts, fragments = frags.37)
s37 <- CreateSeuratObject(s37_assay, assay = "ATAC", meta.data = md.37)

s40_assay <- CreateChromatinAssay(s40.counts, fragments = frags.40)
s40 <- CreateSeuratObject(s40_assay, assay = "ATAC", meta.data = md.40)

s43_assay <- CreateChromatinAssay(s43.counts, fragments = frags.43)
s43 <- CreateSeuratObject(s43_assay, assay = "ATAC", meta.data = md.43)

s45_assay <- CreateChromatinAssay(s45.counts, fragments = frags.45)
s45 <- CreateSeuratObject(s45_assay, assay = "ATAC", meta.data = md.45)

s46_assay <- CreateChromatinAssay(s46.counts, fragments = frags.46)
s46 <- CreateSeuratObject(s46_assay, assay = "ATAC", meta.data = md.46)

s47_assay <- CreateChromatinAssay(s47.counts, fragments = frags.47)
s47 <- CreateSeuratObject(s47_assay, assay = "ATAC", meta.data = md.47)

s50_assay <- CreateChromatinAssay(s50.counts, fragments = frags.50)
s50 <- CreateSeuratObject(s50_assay, assay = "ATAC", meta.data = md.50)

s52_assay <- CreateChromatinAssay(s52.counts, fragments = frags.52)
s52 <- CreateSeuratObject(s52_assay, assay = "ATAC", meta.data = md.52)

s58_assay <- CreateChromatinAssay(s58.counts, fragments = frags.58)
s58 <- CreateSeuratObject(s58_assay, assay = "ATAC", meta.data = md.58)

s66_assay <- CreateChromatinAssay(s66.counts, fragments = frags.66)
s66 <- CreateSeuratObject(s66_assay, assay = "ATAC", meta.data = md.66)

s82_assay <- CreateChromatinAssay(s82.counts, fragments = frags.82)
s82 <- CreateSeuratObject(s82_assay, assay = "ATAC", meta.data = md.82)

s90_assay <- CreateChromatinAssay(s90.counts, fragments = frags.90)
s90 <- CreateSeuratObject(s90_assay, assay = "ATAC", meta.data = md.90)

s96_assay <- CreateChromatinAssay(s96.counts, fragments = frags.96)
s96 <- CreateSeuratObject(s96_assay, assay = "ATAC", meta.data = md.96)




