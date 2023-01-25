

#qc of the object

#mitochondrial
total_counts_per_cell <- colSums(seurat_object@assays$RNA@counts)
mito_genes <- rownames(seurat_object)[grep("^MT-", rownames(seurat_object))]
seurat_object$percent_mito <- colSums(seurat_object@assays$RNA@counts[mito_genes, ])/total_counts_per_cell

#head(mito_genes, 10)

#ribosomal
ribo_genes <- rownames(seurat_object)[grep("^RP[SL]", rownames(seurat_object))]
seurat_object$percent_ribo <- colSums(seurat_object@assays$RNA@counts[ribo_genes, ])/total_counts_per_cell

#hemoglobin
# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
seurat_object <- PercentageFeatureSet(seurat_object, "^HB[^(P)]", col.name = "percent_hb")
seurat_object <- PercentageFeatureSet(seurat_object, "PECAM1|PF4", col.name = "percent_plat")


feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(seurat_object, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
    NoLegend()
    
FeatureScatter(seurat_object, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
