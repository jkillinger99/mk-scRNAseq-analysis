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

# QC and filtering
res4aggr[["percent.mt"]] <- PercentageFeatureSet(res4aggr, pattern = "^mt-")
res4aggr_filtered <- subset(res4aggr, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 25)

# normalization
res4aggr_filtered <- SCTransform(res4aggr_filtered, verbose = FALSE) # could regress out vars using vars.to.regress parameter
res4aggr_filtered <- RunPCA(res4aggr_filtered, verbose = FALSE)
ElbowPlot(res4aggr_filtered, ndims = 50)

res4aggr_filtered <- RunUMAP(res4aggr_filtered, dims = 1:50)
res4aggr_filtered <- FindNeighbors(res4aggr_filtered, dims = 1:50)
res4aggr_filtered <- FindClusters(res4aggr_filtered)

# plotting umap
res4_umap_unintegrated <- DimPlot(res4aggr_filtered, reduction = "umap", group.by = c("seurat_clusters", "Condition", "sample_id"), label = TRUE)


# finding markers by seurat cluster
res4aggr_filtered <- PrepSCTFindMarkers(res4aggr_filtered, verbose = T)

Idents(res4aggr_filtered) = "SCT_snn_res.0.8"
DefaultAssay(res4aggr_filtered) = "SCT"

de_allClusters = FindAllMarkers(res4aggr_filtered, test.use = "wilcox", only.pos = T, min.pct = 0.1)

# small subset of markers by seurat cluster
de_allClusters %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

# heatmap of top markers by seurat cluster
top20<-de_allClusters %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

DoHeatmap(res4aggr_filtered, features = top20$gene) + NoLegend()

hsc=c("Avp", "Klf2", "Crhbp", "Hopx", "Ids", "Bst2", "Spink2")
mep=c("Apoc1","Hbd","Cnr1p1","Klf1","Blvrb","Linc00152","Tmem14c","Itga2b","Myc","Slc40a1")
mkp1=c("Gata2","Fcer1a","Pbx1","Cytl1","Slc40a1")
mkp2=c("Ppbp","Pf4","Nrgn","Plek","Rgs18","Sdrp","Pdlim1")

chat=c("Ppbp","Cxcl7","Cxcl3","Itga2b","Gp9")
chat2=c("Itgb3", "Mpl", "Fli1", "Runx1", "Pf4", "Ppbp", "Treml1")

VlnPlot(res4aggr_filtered, features=chat2)
FeaturePlot(res4aggr_filtered,features=contam)

grep("(?i)^ppbp$|^cxcl7$|^pf4$|^cxcl4$", rownames(res4aggr_filtered), perl=TRUE)
GetAssayData(object = res4aggr_filtered, slot="counts")["Pf4",]


# finding markers by condition
res4aggr_filtered <- PrepSCTFindMarkers(res4aggr_filtered, verbose = T)

Idents(res4aggr_filtered) = "Condition"
DefaultAssay(res4aggr_filtered) = "SCT"

de_Conditions = FindAllMarkers(res4aggr_filtered, test.use = "wilcox", only.pos = T, min.pct = 0.1)

top20C<-de_Conditions %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() 

DoHeatmap(res4aggr_filtered, features = top20C$gene) + NoLegend()
VlnPlot(res4aggr_filtered, features=mkp2)
