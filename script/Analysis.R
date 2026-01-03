############################################
# Melanoma scRNA-seq Analysis with
# Automatic Cluster Annotation (SingleR)
############################################

#-----------------------------
# Load libraries
#-----------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#-----------------------------
# 1. Load and preprocess data
#-----------------------------
gex_data <- Read10X(data.dir = "D:/Melanoma/filtered_feature_bc_matrix/")

seurat_obj <- CreateSeuratObject(
  counts = gex_data,
  project = "Melanoma_Tumor",
  min.cells = 3,
  min.features = 200
)

# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

# QC filtering
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 250 & percent.mt < 10
)

#-----------------------------
# 2. Normalization & feature selection
#-----------------------------
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)

#-----------------------------
# 3. Dimensionality reduction & clustering
#-----------------------------
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj, ndims = 50)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.4) +
  ggtitle("UMAP – Unannotated Clusters")

#-----------------------------
# 4. Automatic cell-type annotation (SingleR)
#-----------------------------

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj)

# Load human reference
ref <- HumanPrimaryCellAtlasData()

# Run SingleR
singler_results <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

# Add SingleR labels to Seurat metadata
seurat_obj$SingleR_label <- singler_results$labels

#-----------------------------
# 5. Automatic cluster renaming
#-----------------------------

# Determine dominant SingleR label per cluster
cluster_annotation <- seurat_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(
    cell_type = names(
      sort(table(SingleR_label), decreasing = TRUE)
    )[1]
  )

# Create named vector: cluster -> cell type
new_cluster_ids <- setNames(
  cluster_annotation$cell_type,
  cluster_annotation$seurat_clusters
)

# Rename cluster identities
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)

# Store final annotation
seurat_obj$cell_type <- Idents(seurat_obj)

# Visualize annotated UMAP
DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) + ggtitle("UMAP – Automatic Cell-Type Annotation (SingleR)")

#-----------------------------
# 6. Differential expression (cell-type based)
#-----------------------------
cluster_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Top 5 markers per cell type
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(5, wt = avg_log2FC)

# Heatmap
DoHeatmap(
  seurat_obj,
  features = top_markers$gene
) + ggtitle("Top Marker Genes per Cell Type")

#-----------------------------
# 7. Functional enrichment analysis
#-----------------------------

# Example: enrichment for one cell type
cell_type_of_interest <- unique(cluster_markers$cluster)[1]

genes_interest <- cluster_markers %>%
  filter(cluster == cell_type_of_interest & p_val_adj < 0.05) %>%
  pull(gene)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(
  genes_interest,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# GO enrichment
ego <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG enrichment
ekegg <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

#-----------------------------
# 8. Visualize enrichment results
#-----------------------------
dotplot(ego, showCategory = 10) +
  ggtitle(paste("GO Enrichment –", cell_type_of_interest))

barplot(ekegg, showCategory = 10) +
  ggtitle(paste("KEGG Pathways –", cell_type_of_interest))

#-----------------------------
# 9. Inspect results
#-----------------------------
head(cluster_markers)
head(ego)
head(ekegg)
