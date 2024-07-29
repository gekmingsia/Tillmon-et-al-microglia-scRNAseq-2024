# loading libraries
library(Seurat) # v5.0.3
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(speckle) # v1.2.0
library(sctransform)
library(glmGamPoi)

# Loading Seurat object
stress.merged.filt.sct.cca.umap <- readRDS("~/Downloads/scRNA/stress_merged_filt_sct_cca_umap.rds")

# Creating group column needed by speckle
stress.merged.filt.sct.cca.umap$group <- as.integer(grepl("cort", stress.merged.filt.sct.cca.umap$orig.ident))

# Checking that column was created
head(stress.merged.filt.sct.cca.umap)
tail(stress.merged.filt.sct.cca.umap)

# Performing differential abundance analysis with Speckle
propeller(clusters=stress.merged.filt.sct.cca.umap$seurat_clusters, sample=stress.merged.filt.sct.cca.umap$orig.ident, group=stress.merged.filt.sct.cca.umap$group)

# Plotting cluster abundances in bar chart
plotCellTypeProps(clusters = stress.merged.filt.sct.cca.umap$seurat_clusters, sample=stress.merged.filt.sct.cca.umap$orig.ident)

# Getting cell proportions for plotting
props <- getTransformedProps(clusters = stress.merged.filt.sct.cca.umap$seurat_clusters, sample = stress.merged.filt.sct.cca.umap$orig.ident, transform = "logit")
head(props)

