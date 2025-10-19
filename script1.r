# load libraries
library(hdf5r)
library(tidyverse)
library(Seurat)
library(patchwork)

setwd("/Users/jackkillinger/Desktop/")

##### res4 data #####
# the non-normalized aggr data

# read data
data <- Read10X_h5("/Users/jackkillinger/Desktop/res4_filtered_feature_bc_matrix.h5")

# format as seurat object
res4aggr <- CreateSeuratObject(counts = data, project = "mk_analysis", min.cells = 3, min.features = 200)
res4aggr[["percent.mt"]] <- PercentageFeatureSet(res4aggr, pattern = "^mt-")

# filter data
res4aggr_filtered <- subset(res4aggr, subset = nFeature_RNA > 500 & nFeature_RNA < 15000 & percent.mt < 25)

# normalize, find variable features, scale, run pca
res4aggr_filtered <- NormalizeData(res4aggr_filtered)
res4aggr_filtered <- FindVariableFeatures(res4aggr_filtered)
res4aggr_filtered <- ScaleData(res4aggr_filtered)
res4aggr_filtered <- RunPCA(res4aggr_filtered)
ElbowPlot(res4aggr_filtered, ndims = 50)

# clustering and umap
res4aggr_filtered <- FindNeighbors(res4aggr_filtered, dims = 1:50)
res4aggr_filtered <- FindClusters(res4aggr_filtered)
res4aggr_filtered <- RunUMAP(res4aggr_filtered, dims = 1:50)

# plot
res4_umap <- DimPlot(res4aggr_filtered, reduction = "umap", label = TRUE)

##### res3 data #####
# the normalized aggr data

# read data
data <- Read10X_h5("/Users/jackkillinger/Desktop/res3_filtered_feature_bc_matrix.h5")

# format as seurat object
res3aggr <- CreateSeuratObject(counts = data, project = "mk_analysis", min.cells = 3, min.features = 200)
res3aggr[["percent.mt"]] <- PercentageFeatureSet(res3aggr, pattern = "^mt-")

# filter data
res3aggr_filtered <- subset(res3aggr, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)

# normalize, find variable features, scale, run pca
res3aggr_filtered <- NormalizeData(res3aggr_filtered)
res3aggr_filtered <- FindVariableFeatures(res3aggr_filtered)
res3aggr_filtered <- ScaleData(res3aggr_filtered)
res3aggr_filtered <- RunPCA(res3aggr_filtered)
ElbowPlot(res3aggr_filtered, ndims = 50)

# clustering and umap
res3aggr_filtered <- FindNeighbors(res3aggr_filtered, dims = 1:50)
res3aggr_filtered <- FindClusters(res3aggr_filtered)
res3aggr_filtered <- RunUMAP(res3aggr_filtered, dims = 1:50)

# plot
res3_umap <- DimPlot(res3aggr_filtered, reduction = "umap", label = TRUE)

res3_umap + res4_umap


##### comparing markers between res3 and res4 #####
# this is from ChatGPT
# ---- 0) Setup ---------------------------------------------------------------
library(Seurat)
suppressPackageStartupMessages({
  if (!requireNamespace("mclust", quietly=TRUE)) install.packages("mclust")
  if (!requireNamespace("aricode", quietly=TRUE)) install.packages("aricode")
})
library(mclust)   # ARI
library(aricode)  # NMI

# Make sure identities exist
Idents(res3aggr_filtered) <- "seurat_clusters"
Idents(res4aggr_filtered) <- "seurat_clusters"

# ---- 1) Label transfer + agreement -----------------------------------------
# Transfer cluster labels from res3 -> res4
features <- SelectIntegrationFeatures(list(res3aggr_filtered, res4aggr_filtered))
anchors  <- FindTransferAnchors(reference = res3aggr_filtered,
                                query     = res4aggr_filtered,
                                dims = 1:50, features = features, normalization.method = "LogNormalize")
pred     <- TransferData(anchorset = anchors, refdata = Idents(res3aggr_filtered), dims = 1:30)
res4aggr_filtered$pred_from_res3 <- pred$predicted.id

# Confusion, simple accuracy, ARI, NMI
tab <- table(predicted=res4aggr_filtered$pred_from_res3,
             res4=res4aggr_filtered$seurat_clusters)
acc <- sum(diag(prop.table(tab)))  # 1-to-1 agreement on diagonal (after greedy label)
# For ARI/NMI we need vectors of same length:
v1 <- as.character(res4aggr_filtered$pred_from_res3)
v2 <- as.character(res4aggr_filtered$seurat_clusters)
ari <- adjustedRandIndex(v1, v2)
nmi <- NMI(v1, v2)

list(confusion=tab, accuracy=acc, ARI=ari, NMI=nmi)
# if ARI and NMI /geq 0.8, clusters are fine and we're okay to assume similarity btwn normal and non-normal
# in this case they are, so we're okay to proceed with just using res4

##### find markers for each cluster in res4 #####
res4_markers <- FindAllMarkers(res4aggr_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
res4_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>%
  print(n = Inf)

res4_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(res4aggr_filtered, features = top10$gene) + NoLegend()

RidgePlot(res4aggr_filtered, features = c("Vwf", "Pf4"))
VlnPlot(res4aggr_filtered, features = c("Vwf", "Pf4"))
FeaturePlot(res4aggr_filtered, features = c("Vwf", "Pf4", "Itga2b"))
DotPlot(res4aggr_filtered, features=c("Vwf","Pf4", "Itga2b"))


DefaultAssay(res4aggr_filtered)   # ensure "RNA" if you used NormalizeData; use "SCT" if SCTransform
rownames(res4aggr_filtered)[grepl("^(Pf4|Cxcl4|Vwf|Itga2b|Gp9|Gp1ba|Gp6|Tubb1|Nbeal2|Selp|Plek|Mpl|Myl9)$",
                                  rownames(res4aggr_filtered), ignore.case=FALSE)]

mk_markers <- c("Itga2b","Tubb1","Pf4","Cxcl4","Ppbp","Gp9","Gp1ba","Gp6","Myl9","Nbeal2","Selp","Plek","Vwf","Mpl")
res4aggr_filtered <- AddModuleScore(res4aggr_filtered, features=list(mk_markers), name="MKscore")
FeaturePlot(res4aggr_filtered, features="MKscore1")
VlnPlot(res4aggr_filtered, features="MKscore1", group.by="seurat_clusters", pt.size=0)
FeaturePlot(res4aggr_filtered, features=c("Vwf","Pf4","Itga2b","Tubb1","Gp9","Gp1ba"), ncol=3)

panel <- list(
  MK = c("Itga2b","Tubb1","Myl9","Pf4","Cxcl4","Gp9","Gp1ba","Gp6","Nbeal2","Selp","Vwf","Mpl","Plek"),
  Ery = c("Hbb-bs","Hbb-bt","Hba-a1","Klf1","Gypa"),
  Myeloid = c("Lyz2","Lst1","S100a8","S100a9","Csf1r","Mpo","Elane"),
  Lymph = c("Cd3e","Cd4","Cd8a","Ms4a1","Cd79a","Nkg7","Klrb1c"),
  HSPC = c("Kit","Ly6a","Procr","Gata2","Runx1")
)
DotPlot(res4aggr_filtered, features=unique(unlist(panel))) + RotatedAxis()

mk_genes <- panel$MK
res4aggr_filtered <- AddModuleScore(res4aggr_filtered, features=list(mk_genes), name="MKscore")
FeaturePlot(res4aggr_filtered, "MKscore1")
VlnPlot(res4aggr_filtered, "MKscore1", group.by="seurat_clusters", pt.size=0)

# fraction of cells expressing ≥1 MK gene (and ≥3 MK genes) per cluster
E <- GetAssayData(res4aggr_filtered, slot="data")  # log1p CPM
mk_idx <- rownames(E) %in% mk_genes
expr_any <- Matrix::colSums(E[mk_idx, , drop=FALSE] > 0) >= 1
expr_3p  <- Matrix::colSums(E[mk_idx, , drop=FALSE] > 0) >= 3
tapply(expr_any, res4aggr_filtered$seurat_clusters, mean)
tapply(expr_3p,  res4aggr_filtered$seurat_clusters, mean)

##### run the analysis on just the first sample #####
