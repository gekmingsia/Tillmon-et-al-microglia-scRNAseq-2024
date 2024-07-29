# loading libraries
library(Seurat) # v5.0.3
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(sctransform)
library(glmGamPoi)

# Load merged and filtered object
stress.merged.filt <- readRDS("~/Downloads/scRNA/stress_merged_filt.rds")

# Normalizing object with scTransform
stress.merged.filt.sct <- SCTransform(stress.merged.filt, verbose = TRUE)

# Finding principal components
stress.merged.filt.sct <- RunPCA(stress.merged.filt.sct, npcs = 50, verbose = TRUE)

# Integrating layers
stress.merged.filt.sct.cca <- IntegrateLayers(object = stress.merged.filt.sct, method = CCAIntegration, normalization.method = "SCT", verbose = TRUE)

# Saving integrated object
saveRDS(stress.merged.filt.sct.cca, "~/Downloads/scRNA/stress_merged_filt_sct_cca.rds")