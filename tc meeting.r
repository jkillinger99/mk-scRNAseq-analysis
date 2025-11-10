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
install.packages("ggrepel")


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
setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/outs/10NOV25")

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

##### Q 1:: removing and retaining upper threshold #####
## (a) removing upper threshold for nFeature_RNA

set.seed(1234)

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

set.seed(1234)

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

set.seed(1234)

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

set.seed(1234)

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

setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/outs/10NOV25")
write.csv(top20noT, "res4 top 20 features no upper threshold.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)
write.csv(top20T, "res4 top 20 features yes upper threshold.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)

##### Q2:: regress out percent.mt #####
# (c) copying and pasting from Q1 without a threshold and adding regression 

set.seed(1234)

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

set.seed(1234)

res4aggr_filtered_reg <- PrepSCTFindMarkers(res4aggr_filtered_reg, verbose = T)
Idents(res4aggr_filtered_reg) = "SCT_snn_res.0.8"
DefaultAssay(res4aggr_filtered_reg) = "SCT"
de_allClusters_reg = FindAllMarkers(res4aggr_filtered_reg, test.use = "wilcox", only.pos = T, min.pct = 0.1)

top20_reg<- de_allClusters_reg %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/outs/10NOV25")
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

##### continuing with analysis based on above findings #####

# differential gene expression by seurat clusters

set.seed(1234)

res4aggr_filtered_noT <- PrepSCTFindMarkers(res4aggr_filtered_noT, verbose = T)
Idents(res4aggr_filtered_noT) = "SCT_snn_res.0.8"
DefaultAssay(res4aggr_filtered_noT) = "SCT"
de_allClustersnoT = FindAllMarkers(res4aggr_filtered_noT, test.use = "wilcox", only.pos = T, min.pct = 0.1)

top20_cluster_noT<- de_allClustersnoT %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/outs/10NOV25")
write.csv(top20_cluster_noT, "res4 top 20 features no upper threshold.csv", quote=FALSE, row.names = FALSE, col.names=TRUE)

# differential gene expression by condition

set.seed(1234)

res4aggr_filtered_noT <- PrepSCTFindMarkers(res4aggr_filtered_noT, verbose = T)
Idents(res4aggr_filtered_noT) = "Condition"
DefaultAssay(res4aggr_filtered_noT) = "SCT"
de_allConditionsnoT = FindAllMarkers(res4aggr_filtered_noT, test.use = "wilcox", only.pos = T, min.pct = 0.1)

top20_condition_noT<- de_allConditionsnoT %>%
    group_by(Condition) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

# differential gene expression by sample id

set.seed(1234)

res4aggr_filtered_noT <- PrepSCTFindMarkers(res4aggr_filtered_noT, verbose = T)
Idents(res4aggr_filtered_noT) = "sample_id"
DefaultAssay(res4aggr_filtered_noT) = "SCT"
de_allSamplesnoT = FindAllMarkers(res4aggr_filtered_noT, test.use = "wilcox", only.pos = T, min.pct = 0.1)

top20_sample_noT<- de_allSamplesnoT %>%
    group_by(sample_id) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

# use this function to search for specific gene names
genenames <- de_allClustersnoT$gene
# head(genenames)
grep("Pdgfb",genenames, ignore.case = TRUE, value = TRUE)

# MK markers given and used by Tianmeng previously
MK_tc <- c("Cd34", "Itga2b", "Mpl", "Gp5", " Stxbp3", "Plek", "Cd9", "Itgb3", "Vwf", "Pf4", "Fermt3", "Mmrn1", "Selp", "Pdgfb")

# expression of markers across clusters
DotPlot(res4aggr_filtered_noT, features = MK_tc, group.by = "seurat_clusters", cols = c("blue", "red"), dot.scale = 8) + theme(legend.position = "right")

# cell counts across cluster by sample id
table(res4aggr_filtered_noT$seurat_clusters, res4aggr_filtered_noT$sample_id)
table(res4aggr_filtered_noT$seurat_clusters, res4aggr_filtered_noT$Condition)

# expression of markers within the same cell type across conditions

set.seed(1234)

res4aggr_filtered_noT$celltype.cond <- paste(res4aggr_filtered_noT$seurat_clusters, res4aggr_filtered_noT$Condition, sep = "_")
Idents(res4aggr_filtered_noT) <- "celltype.cond"
de.10 <- FindMarkers(res4aggr_filtered_noT, ident.1 = "10_poly", ident.2 = "10_oob", test.use = "wilcox", only.pos = F, min.pct = 0.1)
de.12 <- FindMarkers(res4aggr_filtered_noT, ident.1 = "12_poly", ident.2 = "12_oob", test.use = "wilcox", only.pos = F, min.pct = 0.1)

de.10_expression <- de.10 %>% 
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    ungroup()

de.12_expression <- de.12 %>% 
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    ungroup()

setwd("C:/Users/JRK184/OneDrive - University of Pittsburgh/Shea Lab/Studies/Megakaryocytes/analysis/scRNAseq/outs/10NOV25")
write.csv(de.10_expression, "res4 features by cluster 10.csv", quote=FALSE, row.names = TRUE, col.names=TRUE)
write.csv(de.12_expression, "res4 features by cluster 12.csv", quote=FALSE, row.names = TRUE, col.names=TRUE)

# plots of DEGs
# cluster 10
de.10 <- de.10 %>%
  mutate(sig = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
    p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
    TRUE ~ "NotSig"
  ))

top20_up <- de.10 %>%
  filter(sig == "Up") %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20)

top20_down <- de.10 %>%
  filter(sig == "Down") %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 20)

top_labels <- bind_rows(top20_up, top20_down)

ggplot(de.10, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "gray70")) +
  geom_text_repel(
    data = top_labels,
    aes(label = rownames(top_labels)),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Cluster 10 Poly vs OOB",
    x = "Average Log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Direction"
  ) +
  theme_bw()

# cluster 12
de.12 <- de.12 %>%
  mutate(sig = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
    p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
    TRUE ~ "NotSig"
  ))

top20_up <- de.12 %>%
  filter(sig == "Up") %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20)

top20_down <- de.12 %>%
  filter(sig == "Down") %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 20)

top_labels <- bind_rows(top20_up, top20_down)

ggplot(de.12, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "gray70")) +
  geom_text_repel(
    data = top_labels,
    aes(label = rownames(top_labels)),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Cluster 12 Poly vs OOB",
    x = "Average Log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Direction"
  ) +
  theme_bw()


# use this function to search for specific gene names
genenames <- de_allClustersnoT$gene
# head(genenames)
grep("Gpvi",genenames, ignore.case = TRUE, value = TRUE)

# trying with a different set of markers
MK_tc_mod2 <- c("Cd34", "Itga2b", "Mpl", "Gp5", " Stxbp3", "Plek", "Cd9", "Itgb3", "Vwf", "Pf4", "Fermt3", "Mmrn1", "Selp", "Pdgfb", "Gata1", "Zfpm1")

# expression of markers across clusters
DotPlot(res4aggr_filtered_noT, features = MK_tc_mod2, group.by = "seurat_clusters", cols = c("blue", "red"), dot.scale = 8) + theme(legend.position = "right")

# cell counts across cluster by sample id
table(res4aggr_filtered_noT$seurat_clusters, res4aggr_filtered_noT$sample_id)
table(res4aggr_filtered_noT$seurat_clusters, res4aggr_filtered_noT$Condition)