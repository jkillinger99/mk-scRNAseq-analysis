# load libraries
# install.packages("hdf5r")
# install.packages("Seurat")
options(vsc.plot = TRUE)
library(hdf5r)
library(tidyverse)
library(Seurat)
library(patchwork)

setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/R analyses/19OCT25")

##### goal 1: compare res3 (normalized) vs res4 (non-normalized) aggregate data ##### 
res3 <- Read10X_h5("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/Data analysis/Results 3/filtered_feature_bc_matrix.h5")
res4 <- Read10X_h5("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/Data analysis/Results 4/filtered_feature_bc_matrix.h5")

res3aggr <- CreateSeuratObject(counts = res3, project = "mk_analysis", min.cells = 3, min.features = 200)
res4aggr <- CreateSeuratObject(counts = res4, project = "mk_analysis", min.cells = 3, min.features = 200)

res3aggr[["percent.mt"]] <- PercentageFeatureSet(res3aggr, pattern = "^mt-")
res4aggr[["percent.mt"]] <- PercentageFeatureSet(res4aggr, pattern = "^mt-")

# metadata plots by res3 versus res4
res3aggr_metadata <- VlnPlot(res3aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
res4aggr_metadata <- VlnPlot(res4aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

metaaggr_3and4 <- res3aggr_metadata / res4aggr_metadata
png("metaaggr_3and4.png", width = 800, height = 400)
print(metaaggr_3and4)
dev.off()

# features filtered based on TC's papers
res3aggr_filtered <- subset(res3aggr, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)
res4aggr_filtered <- subset(res4aggr, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)


# normalizing, find variable features, scale, run pca
res3aggr_filtered <- NormalizeData(res3aggr_filtered)
res3aggr_filtered <- FindVariableFeatures(res3aggr_filtered)
res3aggr_filtered <- ScaleData(res3aggr_filtered)
res3aggr_filtered <- RunPCA(res3aggr_filtered)
res3ElbowPlot <-ElbowPlot(res3aggr_filtered, ndims = 50)


res4aggr_filtered <- NormalizeData(res4aggr_filtered)
res4aggr_filtered <- FindVariableFeatures(res4aggr_filtered)
res4aggr_filtered <- ScaleData(res4aggr_filtered)
res4aggr_filtered <- RunPCA(res4aggr_filtered)
res4ElbowPlot <- ElbowPlot(res4aggr_filtered, ndims = 50)

EPs_3and4 <- res3ElbowPlot + res4ElbowPlot
png("res3_res4_elbowplots.png", width = 800, height = 400)
print(EPs_3and4)
dev.off()

# clustering and umap
res3aggr_filtered <- FindNeighbors(res3aggr_filtered, dims = 1:50)
res3aggr_filtered <- FindClusters(res3aggr_filtered)
res3aggr_filtered <- RunUMAP(res3aggr_filtered, dims = 1:50)

res4aggr_filtered <- FindNeighbors(res4aggr_filtered, dims = 1:50)
res4aggr_filtered <- FindClusters(res4aggr_filtered)
res4aggr_filtered <- RunUMAP(res4aggr_filtered, dims = 1:50)

res3_umap <- DimPlot(res3aggr_filtered, reduction = "umap", label = TRUE)
res4_umap <- DimPlot(res4aggr_filtered, reduction = "umap", label = TRUE)

umaps_3and4 <- res3_umap + res4_umap
png("res3_res4_umaps.png", width = 800, height = 400)
print(umaps_3and4)
dev.off()

# make cell names unique and add a dataset tag
res3aggr_filtered <- RenameCells(res3aggr_filtered, add.cell.id = "res3")
res4aggr_filtered <- RenameCells(res4aggr_filtered, add.cell.id = "res4")
res3aggr_filtered$dataset <- "res3"
res4aggr_filtered$dataset <- "res4"


features <- SelectIntegrationFeatures(list(res3aggr_filtered, res4aggr_filtered), nfeatures = 3000)
res3aggr_filtered <- ScaleData(res3aggr_filtered, features = features, verbose = FALSE)
res4aggr_filtered <- ScaleData(res4aggr_filtered, features = features, verbose = FALSE)
anchors <- FindIntegrationAnchors(list(res3aggr_filtered, res4aggr_filtered),
                                  anchor.features = features, dims = 1:50)
integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

# reductions
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:50, verbose = FALSE)

# color by res3 vs res4
integrated_plot <- DimPlot(integrated, reduction = "umap", group.by = "dataset", label = TRUE)
png("integrated_plot.png", width = 800, height = 400)
print(integrated_plot)
dev.off()

Idents(integrated) <- "seurat_clusters"
table(integrated$dataset, integrated$seurat_clusters)


##### subgoal: additing groupings to metadata #####
aggr <- read.csv("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/Data analysis/Results 4/aggregation.csv", check.names = FALSE)

aggr$library_id <- as.character(seq_len(nrow(aggr)))
lib_idx <- sub(".*-(\\d+)$", "\\1", colnames(res4aggr))
sample_map <- setNames(aggr$sample_id, aggr$library_id)
cond_map   <- setNames(aggr$Condition,  aggr$library_id)

meta <- data.frame(
  library_id = lib_idx,
  sample_id  = unname(sample_map[lib_idx]),
  Condition  = unname(cond_map[lib_idx]),
  row.names  = colnames(res4aggr),
  stringsAsFactors = FALSE
)

res4aggr <- AddMetaData(res4aggr, metadata = meta)

head(res4aggr@meta.data)

##### goal 2: determine MK-specific gene expression in res4 data #####

# re-running and condensing code above for clarity
res4aggr[["percent.mt"]] <- PercentageFeatureSet(res4aggr, pattern = "^mt-")

# glimpsing at metadata by condition and sample

# retroactively moving these to a "not used" folder and will use density plots instead
res4_cond_nCount_RNA <- VlnPlot(res4aggr, features = "nCount_RNA", layer="counts", group.by="Condition", raster=FALSE, alpha=0.2)
png("/Not used/res4_cond_nCount_RNA.png", width = 800, height = 400)
print(res4_cond_nCount_RNA)
dev.off()
res4_cond_nFeature_RNA <- VlnPlot(res4aggr, features = "nFeature_RNA", layer="counts", group.by="Condition",raster=FALSE,alpha=0.2)
png("/Not used/res4_cond_nFeature_RNA.png", width = 800, height = 400)
print(res4_cond_nFeature_RNA)
dev.off()
res4_cond_percent_mt <- VlnPlot(res4aggr, features = "percent.mt", layer="counts", group.by="Condition",raster=FALSE,alpha=0.2)
png("/Not used/res4_cond_percent_mt.png", width = 800, height = 400)
print(res4_cond_percent_mt)
dev.off()

res4_sample_nCount_RNA <- VlnPlot(res4aggr, features = "nCount_RNA", layer="counts", group.by="sample_id", raster=FALSE, alpha=0.2)
png("/Not used/res4_sample_nCount_RNA.png", width = 800, height = 400)
print(res4_sample_nCount_RNA)
dev.off()
res4_sample_nFeature_RNA <- VlnPlot(res4aggr, features = "nFeature_RNA", layer="counts", group.by="sample_id", raster=FALSE, alpha=0.2)
png("/Not used/res4_sample_nFeature_RNA.png", width = 800, height = 400)
print(res4_sample_nFeature_RNA)
dev.off()
res4_sample_percent_mt <- VlnPlot(res4aggr, features = "percent.mt", layer="counts", group.by="sample_id",raster=FALSE,alpha=0.2)
png("/Not used/res4_sample_percent_mt.png", width = 800, height = 400)
print(res4_sample_percent_mt)
dev.off()

metadata <- res4aggr@meta.data

res4_condition_nCount_RNA_DP <- metadata %>% 
    ggplot(aes(color=Condition, x=nCount_RNA, fill= Condition)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10()
png("res4_condition_nCount_RNA_DP.png", width = 800, height = 400)
print(res4_condition_nCount_RNA_DP)
dev.off()

res4_condition_nFeature_RNA_DP <- metadata %>% 
    ggplot(aes(color=Condition, x=nFeature_RNA, fill= Condition)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 500,color="red",linetype="dotted") +
    geom_vline(xintercept = 5000,color="red",linetype="dotted")
png("res4_condition_nFeature_RNA_DP.png", width = 800, height = 400)
print(res4_condition_nFeature_RNA_DP)
dev.off()

res4_condition_percent.mt_DP <- metadata %>% 
    ggplot(aes(color=Condition, x=percent.mt, fill= Condition)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 25,color="red",linetype="dotted")
png("res4_condition_percent.mt_DP.png", width = 800, height = 400)
print(res4_condition_percent.mt_DP)
dev.off()

res4_sample_id_nCount_RNA_DP <- metadata %>% 
    ggplot(aes(color=sample_id, x=nCount_RNA, fill= sample_id)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10()
png("res4_sample_id_nCount_RNA_DP.png", width = 800, height = 400)
print(res4_sample_id_nCount_RNA_DP)
dev.off()

res4_sample_id_nFeature_RNA_DP <- metadata %>% 
    ggplot(aes(color=sample_id, x=nFeature_RNA, fill= sample_id)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 500,color="red",linetype="dotted")
png("res4_sample_id_nFeature_RNA_DP.png", width = 800, height = 400)
print(res4_sample_id_nFeature_RNA_DP)
dev.off()

res4_sample_id_percent.mt_DP <- metadata %>% 
    ggplot(aes(color=sample_id, x=percent.mt, fill= sample_id)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 25,color="red",linetype="dotted")
png("res4_sample_id_percent.mt_DP.png", width = 800, height = 400)
print(res4_sample_id_percent.mt_DP)
dev.off()

# filtering, normalizing, finding variable features, scaling, pca, clustering, umap from above, repeated
res4aggr_filtered <- subset(res4aggr, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 25)
res4aggr_filtered <- NormalizeData(res4aggr_filtered)
res4aggr_filtered <- FindVariableFeatures(res4aggr_filtered)
res4aggr_filtered <- ScaleData(res4aggr_filtered)
res4aggr_filtered <- RunPCA(res4aggr_filtered)
res4aggr_filtered <- FindNeighbors(res4aggr_filtered, dims = 1:50)
res4aggr_filtered <- FindClusters(res4aggr_filtered)
res4aggr_filtered <- RunUMAP(res4aggr_filtered, dims = 1:50)
res4_umap <- DimPlot(res4aggr_filtered, reduction = "umap", group.by = c("seurat_clusters", "Condition", "sample_id"), label = TRUE)

png("res4_umap.png", width = 800, height = 400)
print(res4_umap)
dev.off()

# find markers that define seurat clusters
res4_markers <- FindAllMarkers(res4aggr_filtered)

# find markers between conditions
res4_markers_conditions <- FindMarkers(res4aggr_filtered, group.by = "Condition")



