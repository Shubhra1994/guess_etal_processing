library(Seurat)
library(tidyverse)
library(Matrix)
library(scDblFinder)
library(scater)
library(hdf5r)
library(loomR)
library(patchwork)
library(scran)
library(scran)

############convert to single cell experiment############
seurat_data <- Read10X(data.dir = "Patient_17_MDS")                      
seurat_obj_p17_MDS <- CreateSeuratObject(counts = seurat_data$'Gene Expression', min.features = 100, project = "Patient_17_MDS")
seurat_obj_p17_sce_MDS <- as.SingleCellExperiment(seurat_obj_p17_MDS)
saveRDS(seurat_obj_p17_sce_MDS,"seurat_obj_p17_sce_MDS.rds")


############run clustering############
set.seed(101000110)
#clusters <- quickCluster(seurat_obj_p17_sce_MDS)
#clusters <- quickCluster(seurat_obj_p17_sce_MDS,method="hclust")
clusters <- quickCluster(seurat_obj_p17_sce_MDS,method="igraph")

seurat_obj_p17_sce_MDS <- computeSumFactors(seurat_obj_p17_sce_MDS, clusters=clusters)
seurat_obj_p17_sce_MDS <- logNormCounts(seurat_obj_p17_sce_MDS)

##################Variance modelling http://bioconductor.org/books/3.15/OSCA.advanced/doublet-detection.html
dec.mam <- modelGeneVarByPoisson(seurat_obj_p17_sce_MDS)
#top.mam <- getTopHVGs(dec.mam, prop=0.1)
top.mam <- getTopHVGs(dec.mam, prop=0.5)

#pca based on denoised data:
sce.mam <- denoisePCA(seurat_obj_p17_sce_MDS, technical=dec.mam, subset.row=top.mam) #Denoise log-expression data by removing principal components corresponding to technical noise. 
sce.mam <- runTSNE(sce.mam, dimred="PCA")



manno.seurat <- as.Seurat(sce.mam, counts = "counts", data = "logcounts")

##################transfer labe;
bm.query <- manno.seurat  
Healthy <- readRDS(file = "Healthy.rds")


bm.anchors <- FindTransferAnchors(reference = Healthy, query = bm.query, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = bm.anchors, refdata = Healthy$Cell_Type4, dims = 1:30)
bm.query <- AddMetaData(bm.query, metadata = predictions)


Healthy <- RunUMAP(Healthy, dims = 1:30, reduction = "pca", return.model = TRUE)
bm.query <- MapQuery(anchorset = bm.anchors, reference = Healthy, query = bm.query, refdata = list(celltype = "Cell_Type4"), reference.reduction = "pca", reduction.model = "umap")

    



############sctransform############
seurat_obj_p17_MDS

#calculate parameters 
seurat_obj_p17_MDS$log10GenesPerUMI <- log10(seurat_obj_p17_MDS$nFeature_RNA) / log10(seurat_obj_p17_MDS$nCount_RNA)
##Add Mitochondrial Gene percentage to metadata slot
seurat_obj_p17_MDS$mitoRatio <- PercentageFeatureSet(object = seurat_obj_p17_MDS, pattern = "^MT-")

#filter out based on above
seurat_obj_p17_MDS <- subset(x = seurat_obj_p17_MDS,
                         subset= (nCount_RNA >= 500) &
                           (nFeature_RNA >= 250) &
                           (log10GenesPerUMI >= 0.80) &
                           (mitoRatio <= 15))	



#transform
seurat_obj_p17_MDS_sctransform <- SCTransform(seurat_obj_p17_MDS, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.7, verbose = FALSE)

saveRDS(ctrl, "ctrl.rds")

##annotate the sctransformed data
pbmc<- seurat_obj_p17_MDS_sctransform
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
output2 <- pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
write.table(output2, file = "Markers.txt")

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

options(bitmapType='cairo')
png("Feature_Lym_umap.png")
FeaturePlot(pbmc, features = c("CD24", "CD79A", "CD79B", "VPREB1", "DNTT", "IGHM", "GNLY"))
dev.off()

options(bitmapType='cairo')
png("Feature_Mye_umap.png")
FeaturePlot(pbmc, features = c("CD14", "S100A8", "S100A9", "CSTA", "LYZ", "MPO", "AZU"))
dev.off()

options(bitmapType='cairo')
png("Feature_Mye2_umap.png")
FeaturePlot(pbmc, features = c("CLC", "PRG2", "ITGA2B", "CA1", "HBB", "PROM1"))
dev.off()


output3 <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
write.table(output3, file = "output3.txt")
options(bitmapType='cairo')
png("Heatmap.png", width = 800, height = 2400)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()

new.cluster.ids <- c("0_monoblast", "1_APC", "2_monocyte", "3_lymphoid", "4_NKspink", "5_Tcell", "6_earlypremylo",
    "7_tcelllike", "8_TandNK", "9_TandNK2", "10_TandNk2", "11_unknown", "12_pDC","13_Bcell")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
options(bitmapType='cairo')
png("Labelled_umap.png", width = 800, height = 800)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
