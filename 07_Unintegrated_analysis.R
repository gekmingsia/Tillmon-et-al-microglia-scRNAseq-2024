# Loading libraries
library(Seurat) # v5.0.3
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(speckle)
library(sctransform)
library(glmGamPoi)

stress.merged.filt <- readRDS("~/Downloads/scRNA/stress_merged_filt.rds")
stress.merged.filt.sct <- SCTransform(stress.merged.filt, verbose = TRUE)
stress.merged.filt.sct <- RunPCA(stress.merged.filt.sct, npcs = 50, verbose = TRUE)
stress.merged.filt.sct <- PrepSCTFindMarkers(stress.merged.filt.sct)
stress.merged.filt.sct <- FindNeighbors(stress.merged.filt.sct, reduction = "pca", dims = 1:20)
stress.merged.filt.sct <- FindClusters(stress.merged.filt.sct, algorithm = 1, resolution = 0.05, verbose = TRUE)
markers <- FindAllMarkers(stress.merged.filt.sct, assay = "SCT", max.cells.per.ident = 2000)
markers['pct.diff'] <- markers['pct.1'] - markers['pct.2']
markers_filt <- filter(markers, avg_log2FC > 0.70 & p_val_adj < 1e-70)
markers_filt %>% group_by(cluster) %>% top_n(4,pct.diff) -> top4_diff
View(top4_diff)

cluster.ids <- c("Fcrls", "Ddit4", "Taco1", "Apoe", "Top2a")
names(cluster.ids) <- levels(stress.merged.filt.sct)
stress.merged.filt.sct <- RenameIdents(stress.merged.filt.sct, cluster.ids)

stress.merged.filt.sct <- RunUMAP(stress.merged.filt.sct, dims = 1:15)
DimPlot(stress.merged.filt.sct, reduction = "umap")


# Pseudobulk DEG analysis
stress.merged.filt.sct$group <- as.integer(grepl("cort", stress.merged.filt.sct$orig.ident))
propeller(clusters=stress.merged.filt.sct$seurat_clusters, sample=stress.merged.filt.sct$orig.ident, group=stress.merged.filt.sct$group)
plotCellTypeProps(clusters = stress.merged.filt.sct$seurat_clusters, sample=stress.merged.filt.sct$orig.ident)
props <- getTransformedProps(clusters = stress.merged.filt.sct$seurat_clusters, sample = stress.merged.filt.sct$orig.ident, transform = "logit")
head(props)
