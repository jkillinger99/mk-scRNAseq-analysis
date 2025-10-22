#install.packages(c("remotes","Seurat","Matrix","fields","RANN","pROC","cowplot"))
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
#install.packages("DoubletFinder")
#install.packages("sctransform")
#install.packages("tidyverse")
#install.packages("hdf5r")
#install.packages("Seurat")
#install.packages('BiocManager')
#BiocManager::install('glmGamPoi')
# remotes::install_github('immunogenomics/presto')

options(vsc.plot = TRUE)
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
setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/R analyses/20OCT25")
res4_location <- "C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/Data analysis/Results 4/filtered_feature_bc_matrix.h5"
res4_metadatalocation <- "C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/Data analysis/Results 4/aggregation.csv"
setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Lab Notebooks/JK/MKs/R analyses/24OCT25/unintegrated")

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

# QC, filtering, dimensionality reduction, PCA
res4aggr[["percent.mt"]] <- PercentageFeatureSet(res4aggr, pattern = "^mt-")
res4aggr_filtered <- res4aggr  %>% 
    subset(subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 25)  %>% 
    SCTransform(verbose = FALSE, vars.to.regress = "percent.mt")  %>% # changed this to regress out any potential mtdna
    RunPCA(verbose = FALSE)  %>% 
    RunUMAP(dims = 1:50)  %>% 
    FindNeighbors(dims = 1:50)  %>% 
    FindClusters()

# plots for most VF, elbow, pca, umap
top10_VF <- head(VariableFeatures(res4aggr_filtered), 10)
top15_VF <- head(VariableFeatures(res4aggr_filtered), 15)

# export most VF to csv
write.table(top15_VF, "res4 top15 most variable transcripts.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)

VF_unlabeled <- VariableFeaturePlot(res4aggr_filtered) + theme(legend.position="top")
VF_labeled <- LabelPoints(plot = VF_unlabeled, points = top15_VF, repel = TRUE) +
    theme(legend.position = "none")

res4_VF <- VF_unlabeled + VF_labeled

res4_dimloadings <- VizDimLoadings(res4aggr_filtered, nfeatures = 10, dims = 1:5, reduction = "pca")
res4_dimloadings_table <- head(res4aggr_filtered[["pca"]], nfeatures = 10, dims = 1:5)

# export loadings to csv 
write.table(res4_dimloadings_table, "res4 1 to 5 loadings with 10 features from pca.csv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ",")

res4_elbowplot_unintegrated <- ElbowPlot(res4aggr_filtered, ndims = 50)
res4_pca_unintegrated <- DimPlot(res4aggr_filtered, reduction = "pca", 
    group.by = c("seurat_clusters", "Condition", "sample_id"), label = TRUE)
res4_umap_unintegrated <- DimPlot(res4aggr_filtered, reduction = "umap", 
    group.by = c("seurat_clusters", "Condition", "sample_id"), label = TRUE)

images <- list(
    elbow = res4_elbowplot_unintegrated,
    pca = res4_pca_unintegrated,
    umap = res4_umap_unintegrated,
    vf = res4_VF,
    dim = res4_dimloadings
)

# for (name in names(images)) {
#   ggsave(
#     filename = paste0("res4_", name, "_unintegrated.png"),
#     plot = images[[name]],
#     width = 20,   
#     height = 10,  
#     units = "in", 
#     dpi = 300     
#   )
# }

# finding markers by seurat cluster
res4aggr_filtered <- PrepSCTFindMarkers(res4aggr_filtered, verbose = T)

Idents(res4aggr_filtered) = "SCT_snn_res.0.8"
DefaultAssay(res4aggr_filtered) = "SCT"

de_allClusters = FindAllMarkers(res4aggr_filtered, test.use = "wilcox", only.pos = T, min.pct = 0.1)

# small subset of markers by seurat cluster
de_allClusters %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

# heatmap of top 20 markers by seurat cluster
top20 <- de_allClusters %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

top5 <- de_allClusters  %>% 
    group_by(cluster)  %>% 
    dplyr::filter(avg_log2FC > 1)  %>% 
    slice_head(n = 5)  %>% 
    ungroup()

top2 <- de_allClusters  %>% 
    group_by(cluster)  %>% 
    dplyr::filter(avg_log2FC > 1)  %>% 
    slice_head(n = 2)  %>% 
    ungroup()

write.table(top20, "res4 top 20 features defining seurat clusters.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)

DoHeatmap(res4aggr_filtered, features = top5$gene, group.by = "seurat_clusters")
DotPlot(res4aggr_filtered, features = top2$gene, group.by = "seurat_clusters") + guides(x = guide_axis(angle = 45)) + coord_flip() + theme_bw()

# finding markers by condition
res4aggr_filtered <- PrepSCTFindMarkers(res4aggr_filtered, verbose = T)

Idents(res4aggr_filtered) = "Condition"
DefaultAssay(res4aggr_filtered) = "SCT"

de_Conditions = FindAllMarkers(res4aggr_filtered, test.use = "wilcox", only.pos = T, min.pct = 0.1)

?FindAllMarkers()

top20C <-de_Conditions %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

top5C <-de_Conditions %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 5) %>%
    ungroup() 

top2C <-de_Conditions %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 2) %>%
    ungroup() 

write.table(top20, "res4 top 20 features defining oob v poly.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)

DoHeatmap(res4aggr_filtered, features = unique(top20C$gene), group.by = "Condition")
DotPlot(res4aggr_filtered, features = unique(top20C$gene), group.by = "Condition") + coord_flip() + theme_bw()

# re-plotting the same answers as above but with some known markers from papers
# as rida suggested, could validate with our own markers

BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scRNAseq")
library(SingleR)
library(celldex)

searchReferences("mouse")

# gene reference library

ref_mouse1 <- fetchReference("mouse_rnaseq","2024-02-26")
ref_mouse2 <- fetchReference("immgen", "2024-02-26")

DefaultAssay(res4aggr_filtered) <- "SCT"             # you already use SCT
sce_test <- GetAssayData(res4aggr_filtered, slot="data")  # log-normalized (SCT) matrix

ref <- MouseRNAseqData()   # or: ref <- ImmGenData() ; try both and compare
pred <- SingleR(test = sce_test, ref = ref, labels = ref$label.main)

res4aggr_filtered$SingleR_label <- pred$labels
DimPlot(res4aggr_filtered, group.by = "SingleR_label", label = TRUE, repel = TRUE)


ref2 <- ImmGenData()   # or: ref <- ImmGenData() ; try both and compare
pred <- SingleR(test = sce_test, ref = ref2, labels = ref2$label.main)

res4aggr_filtered$SingleR_label <- pred$labels
DimPlot(res4aggr_filtered, group.by = "SingleR_label", label = TRUE, repel = TRUE)

# Optional QC: see how SingleR labels align with clusters
table(Cluster = res4aggr_filtered$seurat_clusters, Label = res4aggr_filtered$SingleR_label)

# Optional QC: see how SingleR labels align with clusters
table(Cluster = res4aggr_filtered$seurat_clusters, Label = res4aggr_filtered$SingleR_label)

# searching names of genes

genenames <- de_allClusters$gene
head(genenames)

grep("gypc",genenames, ignore.case = TRUE, value = TRUE)

# some from DOI 10.1182/blood.2021010697.
mk <- c("Vwf", "Pf4", "Itga2b", "Gp9", "Itgb3", "Gata2") 

# from ref (8), and include HSC, MEP, CMP, GMP
precursors <- c("Prss57","Sox4","Hus1","Cdk6", "Egfl7","Cd34","Smim24") # Database does not have: SPINK2 (used "Hus1"), RP11-620J15.3M, KIA0125, CYTL1

# from ref (5)

# rbc

rbc <- c("Hba-a1", "Hba-a2", "Gypc")

# swap out features mapping

FeaturePlot(res4aggr_filtered, features = rbc)

# ref (8)