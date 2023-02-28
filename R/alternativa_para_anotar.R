library(SingleR)
library(dplyr)

so <- combined

so_reference <- readRDS("../morabito_neurons/Anotacion_finalizada.rds")
sce_reference <- as.SingleCellExperiment(so_reference, assay = "RNA")
colData(sce_reference) <- as.data.frame(colData(sce_reference)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce_reference))

sce <- as.SingleCellExperiment(so, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce))

library(scater)

pred <- SingleR(test = sce, ref = sce_reference,
                labels = sce_reference$tags, assay.type.test=1,
                BPPARAM= BiocParallel::MulticoreParam(4)) # 8 CPUs.
pred_modf <- pred[!duplicated(row.names(pred)),]
