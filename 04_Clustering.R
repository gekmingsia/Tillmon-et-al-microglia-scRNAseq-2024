# loading libraries
library(Seurat) # v5.0.3
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(speckle)
library(sctransform)
library(glmGamPoi)

# Loading integrated Seurat object
stress.merged.filt.sct.cca <- readRDS("~/Downloads/scRNA/stress_merged_filt_sct_cca.rds")

# Visualizing PCs
ElbowPlot(stress.merged.filt.sct.cca, ndims = 50, reduction = "pca")
# Choosing 15 PCs for subsequent steps

# Clustering
stress.merged.filt.sct.cca <- FindNeighbors(stress.merged.filt.sct.cca, reduction = "integrated.dr", dims = 1:15)
stress.merged.filt.sct.cca <- FindClusters(stress.merged.filt.sct.cca, algorithm = 1, resolution = 0.07, verbose = TRUE)

# Creating table of marker genes for clusters
stress.merged.filt.sct.cca <- PrepSCTFindMarkers(stress.merged.filt.sct.cca)
markers <- FindAllMarkers(object = stress.merged.filt.sct.cca, assay = "SCT", only.pos = TRUE, max.cells.per.ident = 3000)

# Creating column of difference in percentage of cells expressing marker gene in this cluster vs all other clusters
markers['pct.diff'] <- markers['pct.1'] - markers['pct.2']

# Viewing table of marker genes
View(markers)

# Saving table of marker genes
saveRDS(markers, "~/Downloads/scRNA/markers.rds")

# Creating table of marker genes filtered by expression level and p values
markers_filt <- filter(markers, avg_log2FC > 0.70 & p_val_adj < 1e-70)

# Sorting top 3 marker genes in each cluster based on difference in percentage of cells expressing marker in this vs all other clusters
markers_filt %>% group_by(cluster) %>% top_n(3,pct.diff) -> top3_diff

# Viewing table of sorted genes
View(top3_diff)

# Labeling clusters in Seurat object
cluster.ids <- c("Fcrls", "Apoe", "Taco1", "Ifitm3", "Ccl4")
names(cluster.ids) <- levels(stress.merged.filt.sct.cca)
stress.merged.filt.sct.cca <- RenameIdents(stress.merged.filt.sct.cca, cluster.ids)

# Generating UMAP plot of clusters
stress.merged.filt.sct.cca <- RunUMAP(stress.merged.filt.sct.cca, dims = 1:15, reduction = "integrated.dr")

# Saving Seurat object with UMAP
saveRDS(stress.merged.filt.sct.cca, "~/Downloads/scRNA/stress_merged_filt_sct_cca_umap.rds")

# Plotting UMAP
DimPlot(stress.merged.filt.sct.cca, reduction = "umap")
DimPlot(stress.merged.filt.sct.cca, reduction = "umap", split.by = "orig.ident")

# Plotting heatmap of top 3 marker genes for each cluster
DoHeatmap(stress.merged.filt.sct.cca, features = top3_diff$gene, raster = FALSE)

# Plotting violin plots of marker genes for each cluster
VlnPlot(stress.merged.filt.sct.cca, "Apoe", pt.size = 0)
VlnPlot(stress.merged.filt.sct.cca, "Fcrl2", pt.size = 0)
VlnPlot(stress.merged.filt.sct.cca, "Taco1", pt.size = 0)
VlnPlot(stress.merged.filt.sct.cca, "Ifitm3", pt.size = 0)
VlnPlot(stress.merged.filt.sct.cca, "Ccl4", pt.size = 0)

# Plotting marker genes onto UMAP plot
FeaturePlot(stress.merged.filt.sct.cca, features = "Apoe")
FeaturePlot(stress.merged.filt.sct.cca, features = "Fcrls")
FeaturePlot(stress.merged.filt.sct.cca, features = "Taco1")
FeaturePlot(stress.merged.filt.sct.cca, features = "Ifitm3")
FeaturePlot(stress.merged.filt.sct.cca, features = "Ccl4")

# Plotting microglial markers onto UMAP plot
FeaturePlot(stress.merged.filt.sct.cca, features = "C1qa")
FeaturePlot(stress.merged.filt.sct.cca, features = "Cx3cr1")
FeaturePlot(stress.merged.filt.sct.cca, features = "P2ry12")
FeaturePlot(stress.merged.filt.sct.cca, features = "Tmem119")

# Plotting macrophage markers onto UMAP plot
FeaturePlot(stress.merged.filt.sct.cca, features = "F13a1")
FeaturePlot(stress.merged.filt.sct.cca, features = "Mgl2")
FeaturePlot(stress.merged.filt.sct.cca, features = "Ms4a7")
FeaturePlot(stress.merged.filt.sct.cca, features = "Pf4")

# Plotting monocyte markers onto UMAP plot
FeaturePlot(stress.merged.filt.sct.cca, features = "Ccr2")
FeaturePlot(stress.merged.filt.sct.cca, features = "Ly6c2")

# Creating table of cells in each cluster
cellnumbers <- table(Idents(stress.merged.filt.sct.cca), stress.merged.filt.sct.cca$orig.ident)
View(cellnumbers)

# Generating table of positive and negative marker genes for each cluster
markers <- FindAllMarkers(stress.merged.filt.sct.cca, assay = "SCT", max.cells.per.ident = 2000)

# Filtering marker genes based on Fc and p values
markers_apoe_pos <- filter(markers, cluster == "Apoe" & avg_log2FC > 0.25 & p_val_adj < 1e-40)
View(markers_apoe_pos)
markers_apoe_neg <- filter(markers, cluster == "Apoe" & avg_log2FC < -0.25 & p_val_adj < 1e-40)
View(markers_apoe_neg)
markers_fcrls_pos <- filter(markers, cluster == "Fcrls" & avg_log2FC > 0.25 & p_val_adj < 1e-40)
View(markers_fcrls_pos)
markers_fcrls_neg <- filter(markers, cluster == "Fcrls" & avg_log2FC < -0.25 & p_val_adj < 1e-40)
View(markers_fcrls_neg)
markers_taco1_pos <- filter(markers, cluster == "Taco1" & avg_log2FC > 0.25 & p_val_adj < 1e-40)
View(markers_taco1_pos)
markers_taco1_neg <- filter(markers, cluster == "Taco1" & avg_log2FC < -0.25 & p_val_adj < 1e-40)
View(markers_taco1_neg)
markers_ifitm3_pos <- filter(markers, cluster == "Ifitm3" & avg_log2FC > 0.25 & p_val_adj < 1e-40)
View(markers_ifitm3_pos)
markers_ifitm3_neg <- filter(markers, cluster == "Ifitm3" & avg_log2FC < -0.25 & p_val_adj < 1e-40)
View(markers_ifitm3_neg)
markers_ccl4_pos <- filter(markers, cluster == "Ccl4" & avg_log2FC > 0.25 & p_val_adj < 1e-40)
View(markers_ccl4_pos)
markers_ccl4_neg <- filter(markers, cluster == "Ccl4" & avg_log2FC < -0.25 & p_val_adj < 1e-40)
View(markers_ccl4_neg)
