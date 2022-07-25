---
title: "scRNA-seq, NIH/NIA/IRP"
author: "Bongsoo Park"
date: "12/22/2020"

library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(quanteda)
library(cowplot)
library(edgeR)

###Mouse data - 2 BreastCancer treatments, 2 WD treatments
setwd("~/Dropbox/Desktop/analysis/single-cell-arya")

Tumors_2<-Read10X("./processed_outs/Sample2-count/outs/filtered_feature_bc_matrix/")
Tumors_3<-Read10X("./processed_outs/Sample3-count/outs/filtered_feature_bc_matrix/")
Tumors_4<-Read10X("./processed_outs/Sample4-count/outs/filtered_feature_bc_matrix/")
Tumors_5<-Read10X("./processed_outs/Sample5-count/outs/filtered_feature_bc_matrix/")
Bcell_1<-Read10X("./processed_outs/SampleB1-count/outs/filtered_feature_bc_matrix/")
Bcell_2<-Read10X("./processed_outs/SampleB2-count/outs/filtered_feature_bc_matrix/")
Mono_1<-Read10X("./processed_outs/SampleM1-count/outs/filtered_feature_bc_matrix/")
Mono_2<-Read10X("./processed_outs/SampleM2-count/outs/filtered_feature_bc_matrix/")

Tumors_2 <- CreateSeuratObject(counts = Tumors_2, project = "WD",min.features=1000,min.cells=3)   ###16406 features, 6661 cells
Tumors_3 <- CreateSeuratObject(counts = Tumors_3, project = "WD",min.features=1000,min.cells=3)   ###15476 features, 2339 cells
Tumors_4 <- CreateSeuratObject(counts = Tumors_4, project = "WD",min.features=1000,min.cells=3)   ###16406 features, 6661 cells
Tumors_5 <- CreateSeuratObject(counts = Tumors_5, project = "WD",min.features=1000,min.cells=3)   ###15476 features, 2339 cells
Bcell_1 <- CreateSeuratObject(counts = Bcell_1, project = "WD",min.features=1000,min.cells=3)   ###16406 features, 6661 cells
Bcell_2 <- CreateSeuratObject(counts = Bcell_2, project = "WD",min.features=1000,min.cells=3)   ###15476 features, 2339 cells
Mono_1 <- CreateSeuratObject(counts = Mono_1, project = "WD",min.features=1000,min.cells=3)   ###16406 features, 6661 cells
Mono_2 <- CreateSeuratObject(counts = Mono_2, project = "WD",min.features=1000,min.cells=3)   ###15476 features, 2339 cells

Tumors_2@meta.data[,"Treatment"]<-"Tumors"
Tumors_3@meta.data[,"Treatment"]<-"Tumors"
Tumors_4@meta.data[,"Treatment"]<-"Tumors"
Tumors_5@meta.data[,"Treatment"]<-"Tumors"
Bcell_1@meta.data[,"Treatment"]<-"Bcell"
Bcell_2@meta.data[,"Treatment"]<-"Bcell"
Mono_1@meta.data[,"Treatment"]<-"Mono"
Mono_2@meta.data[,"Treatment"]<-"Mono"

Tumors_2@meta.data[,"Replicate"]<-"1"
Tumors_3@meta.data[,"Replicate"]<-"2"
Tumors_4@meta.data[,"Replicate"]<-"3"
Tumors_5@meta.data[,"Replicate"]<-"4"
Bcell_1@meta.data[,"Replicate"]<-"1"
Bcell_2@meta.data[,"Replicate"]<-"2"
Mono_1@meta.data[,"Replicate"]<-"1"
Mono_2@meta.data[,"Replicate"]<-"2"

Tumors_2@meta.data[,"ID"]<-"Tumors_2"
Tumors_3@meta.data[,"ID"]<-"Tumors_3"
Tumors_4@meta.data[,"ID"]<-"Tumors_4"
Tumors_5@meta.data[,"ID"]<-"Tumors_5"
Bcell_1@meta.data[,"ID"]<-"Bcell_1"
Bcell_2@meta.data[,"ID"]<-"Bcell_2"
Mono_1@meta.data[,"ID"]<-"Mono_1"
Mono_2@meta.data[,"ID"]<-"Mono_2"

WD1<-merge(x=Tumors_2,y=Tumors_4, add.cell.ids=c("Tumors_2","Tumors_4"),project="ARYA")
WD2<-merge(x=Tumors_4,y=Tumors_5, add.cell.ids=c("Tumors_4","Tumors_5"),project="ARYA")
WD3<-merge(x=Bcell_1,y=Mono_1, add.cell.ids=c("Bcell_1","Mono_1"),project="ARYA")
WD4<-merge(x=Mono_1,y=Mono_2, add.cell.ids=c("Mono_1","Mono_2"),project="ARYA")

WD<-merge(x=WD1,y=WD3)
WD<-merge(x=WD,y=WD3)
WD<-merge(x=WD,y=WD4)
###Normalize and find variable features

###Calculate % mitochondrial
WD[["percent.mt"]] <- PercentageFeatureSet(WD, pattern = "^mt-")
summary(WD[["percent.mt"]])  ###Looks pretty good! Low % mito
WD[["percent.rpl"]] <- PercentageFeatureSet(WD, pattern = "^Rpl-")
summary(WD[["percent.rpl"]])  ###Looks pretty good! Low % mito

WD <- SCTransform(WD, vars.to.regress = "percent.mt", verbose = FALSE)
###Dimension reduction
WD <- NormalizeData(WD)
WD <- FindVariableFeatures(WD)
WD <- RunFastMNN(object.list = SplitObject(WD, split.by = "ID"))
WD <- RunUMAP(WD, reduction = "mnn", dims = 1:30)
WD <- FindNeighbors(WD, reduction = "mnn", dims = 1:30)
WD <- FindClusters(WD)
DimPlot(WD, group.by = c("ID"), ncol = 1)


ElbowPlot(WD, ndims = 100) ###First elbow around 25 PCs, not much additional variance around 50 PCs or so
DimHeatmap(WD, dims = c(1:9), cells = 500, balanced = TRUE)

###How do QC parameters look across samples?

Idents(object=WD) <- "ID"

VlnPlot(WD,features="percent.mt",pt.size = 0.00) ###WD2 has a higher % mito
VlnPlot(WD,features="nCount_RNA",pt.size = 0.00) ###Read distributions are similar
VlnPlot(WD,features="nFeature_RNA",pt.size = 0.00) ###WD1 and WD2 have more features - interesting!

# Seurat QC pipeline
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
FeatureScatter(WD, "nCount_RNA", "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0.5)
FeatureScatter(WD, "nCount_RNA", "nFeature_RNA", group.by = "ID", pt.size = 0.5)

# SC transform testing
WD.list<-SplitObject(object = WD, split.by = "ID")

for (i in 1:length(x = WD.list)) {
  #WD.list[[i]] <- NormalizeData(object = WD.list[[i]], verbose = FALSE)
  WD.list[[i]] <- SCTransform(object = WD.list[[i]], 
                                        verbose = FALSE)
}

WD.features <- SelectIntegrationFeatures(object.list = WD.list, nfeatures = 3000)
WD.list <- PrepSCTIntegration(object.list = WD.list, anchor.features = WD.features, 
                                    verbose = FALSE)

WD.anchors <- FindIntegrationAnchors(object.list = WD.list, normalization.method = "SCT", 
                                           anchor.features = WD.features, verbose = FALSE)
WD.breast <- IntegrateData(anchorset = WD.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)
###Dimension reduction
WD.breast <- RunPCA(WD.breast, verbose = FALSE)
WD.breast <- RunUMAP(WD.breast, dims = 1:30)
plots <- DimPlot(WD.breast, group.by = c("tech", "celltype"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3))
### CCA integration approach

WD.list<-SplitObject(object = WD, split.by = "ID")

for (i in 1:length(x = WD.list)) {
  WD.list[[i]] <- NormalizeData(object = WD.list[[i]], verbose = FALSE)
  WD.list[[i]] <- FindVariableFeatures(object = WD.list[[i]], 
                                             selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

reference.list <- WD.list[c("Tumors_2","Tumors_3","Tumors_4","Tumors_5","Bcell_1","Bcell_2","Mono_1","Mono_2")]
WD.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

WD.breast <- IntegrateData(anchorset = WD.anchors, dims = 1:30)

###Switch to the integrated data for downstream analyses
DefaultAssay(object = WD.breast) <- "integrated"

# Run the standard workflow for visualization and clustering
WD.breast <- ScaleData(object = WD.breast, verbose = FALSE)
WD.breast <- RunPCA(object = WD.breast, npcs = 30, verbose = FALSE)
WD.breast <- RunUMAP(object = WD.breast, reduction = "pca", 
                               dims = 1:30)
p1 <- DimPlot(object = WD.breast, reduction = "umap", group.by = "ID")
p2 <- DimPlot(object = WD.breast, reduction = "umap", group.by = "Treatment", 
              label = TRUE, repel = TRUE) + NoLegend()
plot_grid(p1, p2)
p3 <- DimPlot(object = WD.breast, reduction = "umap", group.by = "seurat_clusters", 
              label = TRUE, repel = TRUE) 
plot(p3)
summary(factor(WD.breast$Treatment))

###Dimension reduction
DimHeatmap(WD.breast, dims = c(1:9), cells = 500, balanced = TRUE)

###Graph based clustering
WD.breast <- FindNeighbors(WD.breast, reduction = "pca", dims = 1:30, nn.eps = 0.5)
WD.breast <- FindClusters(WD.breast, resolution = 0.5, n.start = 10)

WD.breast <- RunTSNE(WD.breast, dims = 1:30, nthreads = 4, max_iter = 2000)

p1 <- DimPlot(WD.breast, reduction = "umap", pt.size = 1, group.by="ID",label=TRUE) + ggtitle(label = "ID")
p2 <- DimPlot(WD.breast, reduction = "umap", pt.size = 1, group.by="seurat_clusters",label=TRUE) + ggtitle(label = "Clusters")
plot(p1)
plot(p2)
p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)
plot_grid(p1, p2)
CombinePlots(plots = list(p1,p2),legend="TRUE")

###Csf1r gene
DotPlot(WD.breast,col.min=0,features=c("Csf1r"))+labs(title="Csf1r")


###Csf1r gene
DotPlot(WD.breast,col.min=0,features=c("Pecr","Tmem195","Ptplad2","1810011H11Rik","Fert2","Tlr4","Pon3","Mr1","Arsg","Fcgr1","Camk1","Fgd4","Sqrdl","Csf3r"))+labs(title="Macrophage markers by ImmuneGene consorthium")

FeaturePlot(WD.breast,  features = c("Csf1r","Ccl2","Ms4a1","Cd19"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "tsne")

FeaturePlot(WD.breast,  features = c("Csf3r","Fgd4","Camk1","Fcgr1","Pon3","Tlr4"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "tsne")

FeaturePlot(WD.breast,  features = c("Cd79a","Cd75b","Cd19"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "tsne")

FeaturePlot(WD.breast,  features = c("Pax5"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "tsne")

FeaturePlot(WD.breast,  features = c("Dntt","Kit","Vpreb1","II2ra","Ighc","Igkv"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "tsne")

FeaturePlot(WD.breast,  features = c("Xrcc5","Gm4878","Slco2b1","Gpr77","Gpr160","P2ry13","Tanc2","Sepn1"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "tsne")

###ProB https://www.pnas.org/content/114/20/E3954
DotPlot(WD.breast,col.min=0,features=c("Gr1","Cd11b","Cd16","Cd32","B220","Cd19"))+labs(title="Pre/Pro-B")
DotPlot(WD.breast,col.min=0,features=c("Dntt","Kit","Vpreb1","II2ra","Ighc","Igkv"))+labs(title="Pre/Pro-B")
DotPlot(WD.breast,col.min=0,features=c("Xrcc5","Gm4878","Slco2b1","Gpr77","Gpr160","P2ry13","Tanc2","Sepn1"))+labs(title="Peritoneal Macrophage")
DotPlot(WD.breast,col.min=0,features=c("Mafb","Itga9","Cmklr1","Fez2","Tspan4","Abcc3","Nr1d1","Ptprm","Ctsf","Tfpi"))+labs(title="Lung Macrophage")
DotPlot(WD.breast,col.min=0,features=c("Marco"))+labs(title="Peritoneal+lung")


# Heatmap- General Macrophage markers
FeaturePlot(WD.breast,  features = c("Itgam","Adgre1","Fcgr1","Csf1r","Cebpa","Cebpb"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, B-cell derived macrophage
FeaturePlot(WD.breast,  features = c("Egr1","Id3","Ier2","Clec4n","Ier3","Slc40a1"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request
FeaturePlot(WD.breast,  features = c("Cd79a","Csf1r","Adgre1"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, Cell cycle #1
FeaturePlot(WD,  features = c("Ccna2","Ccnb1","Cdc20","Cenpf","Plk1","Foxm1","Aurka","Id3"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, Cell cycle #1
FeaturePlot(WD,  features = c("Ccna2","Ccnb1","Cdc20","Hmga2"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, Cell cycle #2
FeaturePlot(WD,  features = c("Cenpw","Cenpv","Cenpf","Plk1"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, Cell cycle #3
FeaturePlot(WD,  features = c("Foxm1","Aurka","Id3"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, Cholesterol genes #1
FeaturePlot(WD,  features = c("Hmgcs1","Hmgcr","Fdps","Cyp51","Pmvk","Ldlr","Dhcr24","Fdft1"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, Cholesterol genes #2
FeaturePlot(WD,  features = c("Pmvk","Ldlr","Dhcr24","Fdft1"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, M2 genbes
FeaturePlot(WD,  features = c("Mrc1","Cd163","Cd74","Cd274","Arg1","Il4r","Ciita","Fcgr1a","Tgfb1"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, Cell proliferation markers
FeaturePlot(WD,  features = c("Mki67","Pcna","Ccnd1","Ccne1","Mybl2","Bub1","E2f1","Mcm2","Cdc2"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, Resident macrophage #1
FeaturePlot(WD,  features = c("Sall1", "Sall3", "Zeb1", "Etv1", "Meis3", 
                              "Zfp697", "Zfp788", "Spic", "Nrih3", "Ifr7", "Id1","Nfe2"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

#Chen's request, Resident macrophage #2
FeaturePlot(WD, features = c("Tcf7l2", "Abtb42", "Erg1", "Runx3", "Zfp444", "Ahr", 
                             "Nfatc2", "Irf9", "Nfatc1", "Pparg",  "Atf5"),
            cols = c("lightgrey", "red"),
            min.cutoff = "q9",
            reduction = "umap")

#Chen's request, Resident macrophage #2
FeaturePlot(WD.breast, features = c("Scd2","Mvd","Fdft1","Idi1"),
            cols = c("lightgrey", "red"),
            min.cutoff = "q9",
            reduction = "umap")

WD.cluster.markers <- FindAllMarkers(WD.breast, only.pos = FALSE, logfc.threshold = 0.5)
write.csv(WD.cluster.markers,file="NIH-NIA-cluster-12-22-2020-breast.csv")


