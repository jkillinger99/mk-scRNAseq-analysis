install.packages(c("remotes","Seurat","Matrix","fields","RANN","pROC","cowplot"))
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
install.packages("sctransform")
install.packages("tidyverse")
install.packages("hdf5r")
install.packages("Seurat")
install.packages('BiocManager')
BiocManager::install('glmGamPoi')
install.packages('devtools')
devtools::install_github('immunogenomics/presto')

options(vsc.plot = FALSE)
library(hdf5r)
library(tidyverse)
library(Seurat)
library(patchwork)
library(sctransform)
library(DoubletFinder)
library(BiocManager)
library(glmGamPoi)
library(presto)

##### desktop file dir #####
# setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/R analyses/20OCT25")
# res4_location <- "C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/Data analysis/Results 4/filtered_feature_bc_matrix.h5"
# res4_metadatalocation <- "C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/Data analysis/Results 4/aggregation.csv"
# setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/R analyses/24OCT25/unintegrated")

##### macbook file dir #####
res4_location <- "/Users/jackkillinger/Desktop/res4_filtered_feature_bc_matrix.h5"
res4_metadatalocation <- "/Users/jackkillinger/Desktop/res4_aggregation.csv"
setwd("/Users/jackkillinger/Desktop/mk analysis 09oct25")

# load data
res4 <- Read10X_h5(res4_location)
aggr <- read.csv(res4_metadatalocation, check.names = FALSE)

# create seurat object
res4aggr <- CreateSeuratObject(counts = res4, project = "mk_analysis", min.cells = 3, min.features = 200)

# amend seurat object metadata to include condition and sample id
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
res4aggr[["percent.mt"]] <- PercentageFeatureSet(res4aggr, pattern = "^mt-")

##### Q 1:: removing and retaining upper threshold #####
## (a) removing upper threshold for nFeature_RNA
res4aggr_filtered_noT <- res4aggr  %>% 
    subset(subset = nFeature_RNA > 500 & percent.mt < 25)  %>% 
    SCTransform(verbose = FALSE)  %>%
    RunPCA(verbose = FALSE)  %>% 
    RunUMAP(dims = 1:30, verbose = FALSE)  %>% 
    FindNeighbors(dims = 1:30, verbose = FALSE)  %>% 
    FindClusters(verbose = FALSE)

top15_VF_noT <- head(VariableFeatures(res4aggr_filtered_noT), 15)
VF_unlabeled_noT <- VariableFeaturePlot(res4aggr_filtered_noT) + theme(legend.position="top")
VF_labeled_noT <- LabelPoints(plot = VF_unlabeled_noT, points = top15_VF_noT, repel = TRUE) +
    theme(legend.position = "none")

res4_umap_noT_cluster <- DimPlot(res4aggr_filtered_noT, reduction = "umap", 
    group.by = "seurat_clusters", label = TRUE)
res4_umap_noT_condition <- DimPlot(res4aggr_filtered_noT, reduction = "umap", 
    group.by = "Condition", label = TRUE)
res4_umap_noT_sample <- DimPlot(res4aggr_filtered_noT, reduction = "umap", 
    group.by = "sample_id", label = TRUE)

res4aggr_filtered_noT <- PrepSCTFindMarkers(res4aggr_filtered_noT, verbose = T)
Idents(res4aggr_filtered_noT) = "SCT_snn_res.0.8"
DefaultAssay(res4aggr_filtered_noT) = "SCT"
de_allClustersnoT = FindAllMarkers(res4aggr_filtered_noT, test.use = "wilcox", only.pos = T, min.pct = 0.1)

top20noT<- de_allClustersnoT %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

## (b) retaining current upper threshold for nFeature_RNA
res4aggr_filtered_T <- res4aggr  %>% 
    subset(subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 25)  %>% 
    SCTransform(verbose = FALSE)  %>%
    RunPCA(verbose = FALSE)  %>% 
    RunUMAP(dims = 1:30, verbose = FALSE)  %>% 
    FindNeighbors(dims = 1:30, verbose = FALSE)  %>% 
    FindClusters(verbose = FALSE)

top15_VF_T <- head(VariableFeatures(res4aggr_filtered_T), 15)   
VF_unlabeled_T <- VariableFeaturePlot(res4aggr_filtered_T) + theme(legend.position="top")
VF_labeled_T <- LabelPoints(plot = VF_unlabeled_T, points = top15_VF_T, repel = TRUE) +
    theme(legend.position = "none")

res4_umap_T_cluster <- DimPlot(res4aggr_filtered_T, reduction = "umap", 
    group.by = "seurat_clusters", label = TRUE)
res4_umap_T_condition <- DimPlot(res4aggr_filtered_T, reduction = "umap", 
    group.by = "Condition", label = TRUE)
res4_umap_T_sample <- DimPlot(res4aggr_filtered_T, reduction = "umap", 
    group.by = "sample_id", label = TRUE)

res4aggr_filtered_T <- PrepSCTFindMarkers(res4aggr_filtered_T, verbose = T)
Idents(res4aggr_filtered_T) = "SCT_snn_res.0.8"
DefaultAssay(res4aggr_filtered_T) = "SCT"
de_allClustersT = FindAllMarkers(res4aggr_filtered_T, test.use = "wilcox", only.pos = T, min.pct = 0.1)

top20T <- de_allClustersT %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

## a versus b
quartz()
VF_labeled_noT + VF_labeled_T
res4_umap_noT_cluster + res4_umap_T_cluster
res4_umap_noT_condition + res4_umap_T_condition
res4_umap_noT_sample + res4_umap_T_sample

write.csv(top20noT, "res4 top 20 features no upper threshold.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)
write.csv(top20noT, "res4 top 20 features yes upper threshold.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)

##### Q2:: regress out percent.mt #####
# (c) copying and pasting from Q1 without a threshold and adding regression 
res4aggr_filtered_reg <- res4aggr  %>% 
    subset(subset = nFeature_RNA > 500 & percent.mt < 25)  %>% 
    SCTransform(verbose = FALSE, vars.to.regress = "percent.mt")  %>%
    RunPCA(verbose = FALSE)  %>% 
    RunUMAP(dims = 1:30, verbose = FALSE)  %>% 
    FindNeighbors(dims = 1:30, verbose = FALSE)  %>% 
    FindClusters(verbose = FALSE)

top15_VF_reg <- head(VariableFeatures(res4aggr_filtered_reg), 15)
VF_unlabeled_reg <- VariableFeaturePlot(res4aggr_filtered_reg) + theme(legend.position="top")
VF_labeled_reg <- LabelPoints(plot = VF_unlabeled_reg, points = top15_VF_reg, repel = TRUE) +
    theme(legend.position = "none")

res4_umap_reg_cluster <- DimPlot(res4aggr_filtered_reg, reduction = "umap", 
    group.by = "seurat_clusters", label = TRUE)
res4_umap_reg_condition <- DimPlot(res4aggr_filtered_reg, reduction = "umap", 
    group.by = "Condition", label = TRUE)
res4_umap_reg_sample <- DimPlot(res4aggr_filtered_reg, reduction = "umap", 
    group.by = "sample_id", label = TRUE)

res4aggr_filtered_reg <- PrepSCTFindMarkers(res4aggr_filtered_reg, verbose = T)
Idents(res4aggr_filtered_reg) = "SCT_snn_res.0.8"
DefaultAssay(res4aggr_filtered_reg) = "SCT"
de_allClusters_reg = FindAllMarkers(res4aggr_filtered_reg, test.use = "wilcox", only.pos = T, min.pct = 0.1)

top20_reg<- de_allClusters_reg %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

write.csv(top20_reg, "res4 top 20 features no upper threshold yes regress mt_percent.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)

# a versus c
quartz()
VF_labeled_noT + VF_labeled_reg
quartz()
res4_umap_noT_cluster + res4_umap_reg_cluster
quartz()
res4_umap_noT_condition + res4_umap_reg_condition
quartz()
res4_umap_noT_sample + res4_umap_reg_sample

###### identification of doublets using doublet finder #####
# example code from github readme
sweep.res.list_kidney <- paramSweep(seu_kidney, PCs = 1:10, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

