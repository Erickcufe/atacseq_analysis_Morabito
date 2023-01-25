library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multiprocess", workers = 5)
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


