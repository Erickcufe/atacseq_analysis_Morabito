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

# Add information to identify dataset of origin
s100$sampleID <- "sample-100"
s101$sampleID <- "sample-101"
s17$sampleID <- "sample-17"
s19$sampleID <- "sample-19"
s22$sampleID <- "sample-22"
s27$sampleID <- "sample-27"
s33$sampleID <- "sample-33"
s37$sampleID <- "sample-37"
s40$sampleID <- "sample-40"
s43$sampleID <- "sample-43"
s45$sampleID <- "sample-45"
s46$sampleID <- "sample-46"
s47$sampleID <- "sample-47"
s50$sampleID <- "sample-50"
s52$sampleID <- "sample-52"
s58$sampleID <- "sample-58"
s66$sampleID <- "sample-66"
s82$sampleID <- "sample-82"
s90$sampleID <- "sample-90"
s96$sampleID <- "sample-96"

# Add information to identify phenotype
s100$disease <- "Control"
s101$disease <- "Control"
s17$disease <- "AD"
s19$disease <- "AD"
s22$disease <- "AD"
s27$disease <- "AD"
s33$disease <- "AD"
s37$disease <- "AD"
s40$disease <- "AD"
s43$disease <- "AD"
s45$disease <- "AD"
s46$disease <- "AD"
s47$disease <- "AD"
s50$disease <- "AD"
s52$disease <- "Control"
s58$disease <- "Control"
s66$disease <- "Control"
s82$disease <- "Control"
s90$disease <- "Control"
s96$disease <- "Control"

s100$tangle_stage <- "Stage_2"
s101$tangle_stage <- "Stage_1"
s17$tangle_stage <- "Stage_6"
s19$tangle_stage <- "Stage_6"
s22$tangle_stage <- "Stage_6"
s27$tangle_stage <- "Stage_6"
s33$tangle_stage <- "Stage_5"
s37$tangle_stage <- "Stage_6"
s40$tangle_stage <- "Stage_6"
s43$tangle_stage <- "Stage_6"
s45$tangle_stage <- "Stage_6"
s46$tangle_stage <- "Stage_5"
s47$tangle_stage <- "Stage_5"
s50$tangle_stage <- "Stage_6"
s52$tangle_stage <- "Stage_2"
s58$tangle_stage <- "Stage_2"
s66$tangle_stage <- "Stage_2"
s82$tangle_stage <- "missing"
s90$tangle_stage <- "Stage_1"
s96$tangle_stage <- "Stage_0"

s100$Batch <- 3
s101$Batch <- 2
s17$Batch <- 1
s19$Batch <- 1
s22$Batch <- 2
s27$Batch <- 1
s33$Batch <- 1
s37$Batch <- 1
s40$Batch <- 2
s43$Batch <- 3
s45$Batch <- 3
s46$Batch <- 1
s47$Batch <- 2
s50$Batch <- 1
s52$Batch <- 2
s58$Batch <- 2
s66$Batch <- 2
s82$Batch <- 2
s90$Batch <- 2
s96$Batch <- 3

s100$sex <- "M"
s101$sex <- "F"
s17$sex <- "F"
s19$sex <- "F"
s22$sex <- "M"
s27$sex <- "M"
s33$sex <- "M"
s37$sex <- "F"
s40$sex <- "M"
s43$sex <- "F"
s45$sex <- "F"
s46$sex <- "M"
s47$sex <- "M"
s50$sex <- "F"
s52$sex <- "M"
s58$sex <- "M"
s66$sex <- "F"
s82$sex <- "M"
s90$sex <- "F"
s96$sex <- "M"

s100$APOE <- "missing"
s101$APOE <- "missing"
s17$APOE <- "e33"
s19$APOE <- "e33"
s22$APOE <- "e34"
s27$APOE <- "e23"
s33$APOE <- "e33"
s37$APOE <- "e34"
s40$APOE <- "e34"
s43$APOE <- "missing"
s45$APOE <- "e34"
s46$APOE <- "e23"
s47$APOE <- "e33"
s50$APOE <- "e33"
s52$APOE <- "e33"
s58$APOE <- "e33"
s66$APOE <- "e33"
s82$APOE <- "missing"
s90$APOE <- "missing"
s96$APOE <- "missing"

# Merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = s100,
  y = list(s101, s17, s19, s22, s27, s33, s37, s40, s43, s45, s46, s47, s50,
           s52, s58, s66, s82, s90, s96),
  add.cell.ids = c("sample_100", "sample_101", "sample_17", "sample_19", "sample_22", "sample_27", "sample_33",
                   "sample_37", "sample_40", "sample_43", "sample_45", "sample_46", "sample_47", "sample_50",
                   "sample_52", "sample_58", "sample_66", "sample_82", "sample_90", "sample_96")
)

combined[["ATAC"]]
saveRDS(combined, "pipeline_snATAC.rds")
# Binarize peaks
combined@assays$ATAC@counts@x[combined@assays$ATAC@counts@x > 0] <- 1

################################################################################
# Step 02: Quality Control
################################################################################

combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$nCount_ATAC * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

# for each sample, compute nucleosome banding
# note: download chromosome sizes file from UCSC https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/
chrom.sizes <- read.csv('../Datos_scRNA/morabito_data/snATAC/hg38.chrom.sizes.txt', sep='\t', header=F)
colnames(chrom.sizes) <- c('chr', 'size')
combined$nucleosome_signal <- NA
combined$nucleosome_group <- NA
samples <- unique(combined$sampleID)
for(i in 1:length(samples)){
  print(samples[i])
  temp <- NucleosomeSignal(
    subset(combined, sampleID == samples[i]),
    region=paste0(chrom.sizes[1,][[1]],"-1-",chrom.sizes[1,][[2]])
  )
  temp$nucleosome_group <- ifelse(temp$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
  combined$nucleosome_signal <- ifelse(combined$sampleID == samples[i], temp$nucleosome_signal, combined$nucleosome_signal)
  combined$nucleosome_group <- ifelse(combined$sampleID == samples[i], temp$nucleosome_group, combined$nucleosome_group)
}
combined$pass_qc <- ifelse(combined$peak_region_fragments > 300 & combined$peak_region_fragments < 10000 & combined$pct_reads_in_peaks > 15 & combined$blacklist_ratio < 0.01 & combined$nucleosome_signal < 10, TRUE, FALSE)
combined <- combined[,combined$pass_qc]

saveRDS(combined, "pipeline_snATAC_QualityPartI_Done.rds")

################################################################################
# Step 03: Primary Processing
################################################################################

library(monocle3)

NucSeq.atac <- combined
# NucSeq.atac <- readRDS("mergeSamples_snATAC_geneActivity_PartIII_Done.rds")
peak_matrix <- GetAssayData(NucSeq.atac, slot='counts')
peaks <- data.frame(as.character(rownames(peak_matrix)))
rownames(peaks) <- rownames(peak_matrix)
peaks <- as.data.frame(cbind(peaks,peaks))
colnames(peaks) <- c("GeneSymbol", "gene_short_name")
NucSeq.atac_cds <- new_cell_data_set(
  peak_matrix,
  cell_metadata=NucSeq.atac@meta.data,
  #gene_metadata=genes
)
NucSeq.atac_cds <- detect_genes(NucSeq.atac_cds)
NucSeq.atac_cds <- estimate_size_factors(NucSeq.atac_cds)
NucSeq.atac_cds <- preprocess_cds(NucSeq.atac_cds, method = "LSI")
NucSeq.atac_cds <- align_cds(NucSeq.atac_cds, preprocess_method='LSI', alignment_group = "Batch")
NucSeq.atac_cds <- reduce_dimension(NucSeq.atac_cds, reduction_method = 'UMAP', preprocess_method = "LSI")
NucSeq.atac_cds <- cluster_cells(NucSeq.atac_cds)
umap_coords <- reducedDims(NucSeq.atac_cds)$UMAP

# alternativa
NucSeq.atac <- RunTFIDF(NucSeq.atac)
NucSeq.atac <- FindTopFeatures(NucSeq.atac, min.cutoff = 20)
NucSeq.atac <- RunSVD(
  object = NucSeq.atac,
  assay = 'ATAC',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
NucSeq.atac <- RunUMAP(NucSeq.atac, reduction='lsi', dims=1:30)

# monocle_umap <- NucSeq.atac_cds@reduce_dim_aux[["UMAP"]][rownames(NucSeq.atac_cds@reduce_dim_aux[["UMAP"]]) %in% colnames(NucSeq.atac),]
# colnames(monocle_umap) <- c('UMAP_1', 'UMAP_2')
# all.equal(rownames(monocle_umap), colnames(NucSeq.atac))
NucSeq.atac@reductions$umap@cell.embeddings <- umap_coords
NucSeq.atac$monocle_clusters_umap <- clusters(NucSeq.atac_cds, reduction_method='UMAP')
monocle_aligned <- reducedDims(NucSeq.atac_cds)$Aligned
rownames(monocle_aligned) <- colnames(NucSeq.atac_cds)
colnames(monocle_aligned) <- paste0('LSI_', seq(1:ncol(monocle_aligned)))
all.equal(rownames(monocle_aligned), colnames(NucSeq.atac))
NucSeq.atac@reductions$lsi <- CreateDimReducObject(
  embeddings=monocle_aligned,
  key="LSI_",
  assay="peaks"
)

saveRDS(NucSeq.atac, "pipeline_snATAC_clusterUMAP_PartII_Done.rds")


################################################################################
# Step 04: Construct gene activity matrix
################################################################################
library(Signac)
combined <- readRDS("pipeline_snATAC_clusterUMAP_PartII_Done.rds")
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
# annotations <- annotations[annotations$]
anno <- GenomeInfoDb::getChromInfoFromUCSC("hg38")

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(combined) <- annotations

gene.activities <- GeneActivity(combined)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- Seurat::NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)

saveRDS(combined, "pipeline_snATAC_geneActivity_PartIII_Done.rds")



library(EnsDb.Hsapiens.v86)
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(NucSeq.atac),
  chunk = 10
)
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]
NucSeq.atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
NucSeq.atac <- NormalizeData(
  object = NucSeq.atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(NucSeq.atac$nCount_RNA)
)
DefaultAssay(NucSeq.atac) <- 'RNA'

