#install.packages(c("remotes","Seurat","Matrix","fields","RANN","pROC","cowplot"))
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
#install.packages("DoubletFinder")
#install.packages("sctransform")
#install.packages("tidyverse")
#install.packages("hdf5r")
#install.packages("Seurat")
#install.packages('BiocManager')
#BiocManager::install('glmGamPoi')
options(vsc.plot = TRUE)
library(hdf5r)
library(tidyverse)
library(Seurat)
library(patchwork)
library(sctransform)
library(DoubletFinder)
library(BiocManager)
library(glmGamPoi)

##### desktop file dir #####
setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/R analyses/20OCT25")
res4_location <- "C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/Data analysis/Results 4/filtered_feature_bc_matrix.h5"
res4_metadatalocation <- "C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/Data analysis/Results 4/aggregation.csv"

# now performing on integrated data
res4 <- Read10X_h5(res4_location)
aggr <- read.csv(res4_metadatalocation, check.names = FALSE)
res4aggr <- CreateSeuratObject(counts = res4, project = "mk_analysis", min.cells = 3, min.features = 200)
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

# spliting data by Condition
res4aggr[["percent.mt"]] <- PercentageFeatureSet(res4aggr, pattern = "^mt-")
res4aggr_filtered <- subset(res4aggr, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 25)
res4aggr_filtered[["RNA"]] <- split(res4aggr_filtered[["RNA"]], f = res4aggr_filtered$Condition) 

res4aggr_filtered <- SCTransform(res4aggr_filtered, verbose = FALSE)
res4aggr_filtered <- RunPCA(res4aggr_filtered, verbose = FALSE)
ElbowPlot(res4aggr_filtered, ndims = 50)

res4aggr_filtered <- RunUMAP(res4aggr_filtered, dims = 1:50)
res4aggr_filtered <- FindNeighbors(res4aggr_filtered, dims = 1:50)
res4aggr_filtered <- FindClusters(res4aggr_filtered)

res4aggr_filtered <- IntegrateLayers(res4aggr_filtered, method = CCAIntegration, normalization.method = "SCT", verbose = FALSE)
res4aggr_filtered <- FindNeighbors(res4aggr_filtered, reduction = "integrated.dr", dims = 1:30)
res4aggr_filtered <- FindClusters(res4aggr_filtered, resolution = 0.6)
res4aggr_filtered <- RunUMAP(res4aggr_filtered, dims = 1:30, reduction = "integrated.dr")
DimPlot(res4aggr_filtered, reduction = "umap", group.by = c("seurat_clusters", "Condition", "sample_id"), label = TRUE)

res4aggr_filtered <- PrepSCTFindMarkers(res4aggr_filtered, verbose = FALSE)
