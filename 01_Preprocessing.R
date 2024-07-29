  # loading libraries
library(Seurat) # v5.0.3
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(sctransform)
library(glmGamPoi)


# Loading data
Ctrl1.data <- Read10X_h5("~/Downloads/scRNA/Ctrl1_filtered_feature_bc_matrix.h5")
Ctrl2.data <- Read10X_h5("~/Downloads/scRNA/Ctrl2_filtered_feature_bc_matrix.h5")
Ctrl3.data <- Read10X_h5("~/Downloads/scRNA/Ctrl3_filtered_feature_bc_matrix.h5")
Ctrl4.data <- Read10X_h5("~/Downloads/scRNA/Ctrl4_filtered_feature_bc_matrix.h5")
Cort1.data <- Read10X_h5("~/Downloads/scRNA/Cort1_filtered_feature_bc_matrix.h5")
Cort2.data <- Read10X_h5("~/Downloads/scRNA/Cort2_filtered_feature_bc_matrix.h5")
Cort3.data <- Read10X_h5("~/Downloads/scRNA/Cort3_filtered_feature_bc_matrix.h5")
Cort4.data <- Read10X_h5("~/Downloads/scRNA/Cort4_filtered_feature_bc_matrix.h5")

# Checking data dimensions
dim(Ctrl1.data)
dim(Ctrl2.data)
dim(Ctrl3.data)
dim(Ctrl4.data)
dim(Cort1.data)
dim(Cort2.data)
dim(Cort3.data)
dim(Cort4.data)

# Creating Seurat objects
Ctrl1 <- CreateSeuratObject(Ctrl1.data, project = "ctrl-1", min.cells = 5)
Ctrl2 <- CreateSeuratObject(Ctrl2.data, project = "ctrl-2", min.cells = 5)
Ctrl3 <- CreateSeuratObject(Ctrl3.data, project = "ctrl-3", min.cells = 5)
Ctrl4 <- CreateSeuratObject(Ctrl4.data, project = "ctrl-4", min.cells = 5)
Cort1 <- CreateSeuratObject(Cort1.data, project = "cort-1", min.cells = 5)
Cort2 <- CreateSeuratObject(Cort2.data, project = "cort-2", min.cells = 5)
Cort3 <- CreateSeuratObject(Cort3.data, project = "cort-3", min.cells = 5)
Cort4 <- CreateSeuratObject(Cort4.data, project = "cort-4", min.cells = 5)

# Checking object dimensions
Ctrl1
Ctrl2
Ctrl3
Ctrl4
Cort1
Cort2
Cort3
Cort4

# Checking columns of objects
head(Ctrl1)
head(Ctrl2)
head(Ctrl3)
head(Ctrl4)
head(Cort1)
head(Cort2)
head(Cort3)
head(Cort4)

# Merging objects
stress.raw <- merge(Ctrl1, c(Ctrl2, Cort1, Cort2, Ctrl3, Ctrl4, Cort3, Cort4), add.cell.ids = c("ctrl_1", "ctrl_2", "cort_1", "cort_2", "ctrl_3", "ctrl_4", "cort_3", "cort_4"))

# Checking merged object layers
stress.raw

# Checking merged object cells
table(stress.raw$orig.ident)

# Saving merged object
saveRDS(stress.raw, "~/Downloads/scRNA/stress_merged.rds")


