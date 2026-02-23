---
title: "scRNA-seq Analysis of Breast Cancer Samples (Improved)"
author: "Bongsoo Park"
date: "2026-02-22"
---

# Improved Single-Cell RNA-Seq Analysis Pipeline

# 1. Setup: Load Libraries
# ----------------------------------------------------------------
# It is good practice to load all required libraries at the beginning.

library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(patchwork) # Modern replacement for CombinePlots
library(argparse) # For command-line arguments

# 2. Argument Parsing
# ----------------------------------------------------------------
# Using argparse to make the script portable and remove hardcoded paths.

parser <- ArgumentParser(description = "Process 10x Genomics scRNA-seq data.")
parser$add_argument("--input_dir", help = "Path to the directory containing sample folders (e.g., ./processed_outs/).", required = TRUE)
parser$add_argument("--output_dir", help = "Path to the directory to save results.", default = ".")
parser$add_argument("--seed", type = "integer", default = 42, help = "Random seed for reproducibility.")

args <- parser$parse_args()

set.seed(args$seed)

# 3. Data Loading and Initial Object Creation
# ----------------------------------------------------------------
# A more robust way to load data by defining sample names and paths.

sample_info <- list(
    list(id = "Tumors_2", path = "Sample2-count/outs/filtered_feature_bc_matrix/", treatment = "Tumors", replicate = "1"),
    list(id = "Tumors_3", path = "Sample3-count/outs/filtered_feature_bc_matrix/", treatment = "Tumors", replicate = "2"),
    list(id = "Tumors_4", path = "Sample4-count/outs/filtered_feature_bc_matrix/", treatment = "Tumors", replicate = "3"),
    list(id = "Tumors_5", path = "Sample5-count/outs/filtered_feature_bc_matrix/", treatment = "Tumors", replicate = "4"),
    list(id = "Bcell_1", path = "SampleB1-count/outs/filtered_feature_bc_matrix/", treatment = "Bcell", replicate = "1"),
    list(id = "Bcell_2", path = "SampleB2-count/outs/filtered_feature_bc_matrix/", treatment = "Bcell", replicate = "2"),
    list(id = "Mono_1", path = "SampleM1-count/outs/filtered_feature_bc_matrix/", treatment = "Mono", replicate = "1"),
    list(id = "Mono_2", path = "SampleM2-count/outs/filtered_feature_bc_matrix/", treatment = "Mono", replicate = "2")
)

all_objects <- lapply(sample_info, function(sample) {
    data <- Read10X(file.path(args$input_dir, sample$path))
    obj <- CreateSeuratObject(counts = data, project = "WD", min.features = 1000, min.cells = 3)
    obj@meta.data$Treatment <- sample$treatment
    obj@meta.data$Replicate <- sample$replicate
    obj@meta.data$ID <- sample$id
    return(obj)
})

# 4. Merging and Quality Control
# ----------------------------------------------------------------
# Correctly merge all samples and perform QC.

# Corrected merge logic
wd_merged <- merge(all_objects[[1]], y = all_objects[2:length(all_objects)], 
                   add.cell.ids = sapply(sample_info, `[[`, "id"), 
                   project = "ARYA")

# Calculate mitochondrial and ribosomal protein percentages
wd_merged[["percent.mt"]] <- PercentageFeatureSet(wd_merged, pattern = "^mt-")
wd_merged[["percent.rpl"]] <- PercentageFeatureSet(wd_merged, pattern = "^Rpl-")

# 5. SCTransform Integration (Recommended Method)
# ----------------------------------------------------------------
# Using the modern and recommended SCTransform workflow for integration.

wd_list <- SplitObject(object = wd_merged, split.by = "ID")

wd_list <- lapply(X = wd_list, FUN = SCTransform, vars.to.regress = "percent.mt", verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = wd_list, nfeatures = 3000)
wd_list <- PrepSCTIntegration(object.list = wd_list, anchor.features = features, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = wd_list, normalization.method = "SCT", 
                                  anchor.features = features, verbose = FALSE)

wd_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

# 6. Downstream Analysis (PCA, Clustering, UMAP)
# ----------------------------------------------------------------
# Perform standard downstream analysis on the integrated data.

DefaultAssay(wd_integrated) <- "integrated"

wd_integrated <- RunPCA(wd_integrated, verbose = FALSE)
wd_integrated <- RunUMAP(wd_integrated, dims = 1:30)

wd_integrated <- FindNeighbors(wd_integrated, reduction = "pca", dims = 1:30)
wd_integrated <- FindClusters(wd_integrated, resolution = 0.5)

# 7. Visualization
# ----------------------------------------------------------------
# Generate plots to visualize the results.

# UMAP plots using patchwork
p1 <- DimPlot(wd_integrated, reduction = "umap", group.by = "ID") + ggtitle("UMAP by Sample ID")
p2 <- DimPlot(wd_integrated, reduction = "umap", group.by = "Treatment", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("UMAP by Treatment")
p3 <- DimPlot(wd_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("UMAP by Cluster")

# Combine plots
(p1 | p2) / p3
ggsave(file.path(args$output_dir, "UMAP_summary.png"), width = 12, height = 10)

# Feature plots for key markers
FeaturePlot(wd_integrated, features = c("Csf1r", "Ccl2", "Ms4a1", "Cd19"), 
            cols = c("lightgrey", "red"), min.cutoff = "q9", reduction = "umap")
ggsave(file.path(args$output_dir, "FeaturePlot_Markers.png"), width = 10, height = 8)

# 8. Find Cluster Markers and Save Results
# ----------------------------------------------------------------

all_markers <- FindAllMarkers(wd_integrated, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(all_markers, file = file.path(args$output_dir, "All_Cluster_Markers.csv"))

print("Analysis complete. Results saved to output directory.")

