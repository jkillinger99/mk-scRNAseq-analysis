options(vsc.plot = FALSE)
library(hdf5r)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(patchwork)
library(sctransform)
library(DoubletFinder)
library(BiocManager)
library(glmGamPoi)
library(presto)
library(ggrepel)

##### desktop file dir #####
res4_location <- "C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/results 4/filtered_feature_bc_matrix.h5"
res4_metadatalocation <- "C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/results 4/aggregation.csv"
setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/outs/15NOV25")

##### macbook file dir #####
# res4_location <- "/Users/jackkillinger/Desktop/res4_filtered_feature_bc_matrix.h5"
# res4_metadatalocation <- "/Users/jackkillinger/Desktop/res4_aggregation.csv"
# setwd("/Users/jackkillinger/Desktop/mk analysis 09oct25")

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

set.seed(1234)

# threshold added @ less than 20K features to reduce definitive doublets 

res4aggr_processed <- res4aggr  %>% 
    subset(subset = nFeature_RNA > 500 &  nFeature_RNA < 200000 & percent.mt < 25)  %>% 
    SCTransform(verbose = FALSE)  %>%
    RunPCA(verbose = FALSE)  %>% 
    RunUMAP(dims = 1:30, verbose = FALSE)  %>% 
    FindNeighbors(dims = 1:30, verbose = FALSE)  %>% 
    FindClusters(verbose = FALSE)

res4_umap_T_cluster <- DimPlot(res4aggr_processed, reduction = "umap", 
    group.by = "seurat_clusters", label = TRUE)

set.seed(1234)

# cluster

res4aggr_processed <- PrepSCTFindMarkers(res4aggr_processed, verbose = T)
Idents(res4aggr_processed) = "SCT_snn_res.0.8"
DefaultAssay(res4aggr_processed) = "SCT"
allClusters = FindAllMarkers(res4aggr_processed, test.use = "wilcox", only.pos = F, min.pct = 0.1)

# oob v trauma (Condition)

Idents(res4aggr_processed) = "Condition"
DefaultAssay(res4aggr_processed) = "SCT"
polyVoob = FindAllMarkers(res4aggr_processed, test.use = "wilcox", only.pos = F, min.pct = 0.1)

top20T <- polyVoob %>%
    group_by(Condition) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

# outs

setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/outs/15NOV25")
write.csv(allClusters, "allClusters.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)
write.csv(polyVoob, "polyVoob_nonspecific.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)

# filtering out specific clusters
clusters_to_remove <- c(0, 6, 11, 13, 15, 17)
res4_filt <- subset(res4aggr_processed, subset = !(SCT_snn_res.0.8 %in% clusters_to_remove))

set.seed(1234)

res4_filt <- res4_filt %>%
  SCTransform(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(verbose = FALSE)

DimPlot(res4_filt, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(res4_filt, reduction = "umap", group.by = "Condition", label = TRUE)

# cluster

res4_filt <- PrepSCTFindMarkers(res4_filt, verbose = T)
Idents(res4_filt) = "SCT_snn_res.0.8"
DefaultAssay(res4_filt) = "SCT"
allClusters_filt = FindAllMarkers(res4_filt, test.use = "wilcox", only.pos = F, min.pct = 0.1)

# oob v trauma (Condition)

Idents(res4_filt) = "Condition"
DefaultAssay(res4_filt) = "SCT"
polyVoob_filt = FindAllMarkers(res4_filt, test.use = "wilcox", only.pos = F, min.pct = 0.1)

setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/outs/15NOV25")
write.csv(allClusters_filt, "allClusters_filt.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)
write.csv(polyVoob_filt, "polyVoob_filt.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)