# loading libraries
library(Seurat) # v5.0.3
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(sctransform)
library(glmGamPoi)

# Load merged object after preprocessing
stress.merged <- readRDS("~/Downloads/scRNA/stress_merged.rds")

# Visualize transcripts per cell
VlnPlot(stress.merged, features = "nCount_RNA", group.by = "orig.ident")
VlnPlot(stress.merged, features = "nCount_RNA", group.by = "orig.ident", y.max = 5000)

# Visualize genes per cell
VlnPlot(stress.merged, features = "nFeature_RNA", group.by = "orig.ident")
VlnPlot(stress.merged, features = "nFeature_RNA", group.by = "orig.ident", y.max = 2000)

# Creating and visualizing mitochondrial genes
stress.merged[["percent.mt"]] <- PercentageFeatureSet(stress.merged, pattern = "^mt-")
VlnPlot(stress.merged, features = "percent.mt", group.by = "orig.ident")

# Filtering by mito percentage < 10 and 200 < nFeature_RNA < 7000
stress.merged.filt <- subset(stress.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)

# Saving filtered object
saveRDS(stress.merged.filt, file = "~/Downloads/scRNA/stress_merged_filt.rds")
