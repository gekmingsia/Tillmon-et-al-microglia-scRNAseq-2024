# loading libraries
library(Seurat) # v5.0.3
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(sctransform)
library(glmGamPoi)

# Loading Seurat object
stress.merged.filt.sct.cca.umap <- readRDS("~/Downloads/scRNA/stress_merged_filt_sct_cca_umap.rds")

# Loading Libra libraries
library(limma)
library(edgeR)
library(DESeq2)
library(Libra)

# Creating group column required for Libra, assigning ctrl as group 1, cort as group 0
stress.merged.filt.sct.cca.umap$group <- as.integer(grepl("ctrl", stress.merged.filt.sct.cca.umap$orig.ident))
head(stress.merged.filt.sct.cca.umap)
tail(stress.merged.filt.sct.cca.umap)

# Creating required cell type column for Libra
stress.merged.filt.sct.cca.umap$cell.type <- "microglia"
head(stress.merged.filt.sct.cca.umap)

# Pseudobulk DEGs
degs_sctc <- run_de(stress.merged.filt.sct.cca.umap, meta = "meta.data", cell_type_col = "cell.type", label_col = "group", replicate_col = "orig.ident", de_method = "edgeR")

# Visualizing DEGs with EnhancedVolcano
library(EnhancedVolcano)
EnhancedVolcano(degs_sctc, lab = degs_sctc$gene, x = "avg_logFC", y = "p_val_adj")
EnhancedVolcano(degs_sctc, lab = degs_sctc$gene, x = "avg_logFC", y = "p_val_adj", pCutoff = 10e-4, FCcutoff = 1)

# Visualizing specific genes in volcano plot
EnhancedVolcano(degs_sctc, lab = degs_sctc$gene, x = "avg_logFC", y = "p_val_adj", selectLab = c('Apoe', 'Ifitm3'))

# Gene set enrichment analysis (GSEA) of DEGs with clusterProfiler
library(clusterProfiler)
library(org.Mm.eg.db)
up_de <- degs_sctc %>% filter(avg_logFC > 1, p_val_adj < 10e-4) %>% pull(gene)
up_de_go <- enrichGO(up_de, keyType = "SYMBOL", OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)

# Plotting GSEA results
barplot(up_de_go, showCategory = 10)
dotplot(up_de_go, showCategory = 10, orderBy = "p.adjust", decreasing = FALSE, font.size = 10)


