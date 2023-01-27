####Mostly Copied from https://github.com/hbctraining/scRNA-seq_online/tree/master/lessons

##The 3 (barcode, feature, matrix) scRNA-seq files should be stored in scRNA_dir/(name of scRNA-seq sample)/

library(Seurat)
library(tidyverse)
library(Matrix)

scRNA_dir <- "/home/msakuma/scRNA_seq/"
scRNA_file <- c("Patient_17_MDS", "Patient_17_sAML", "Patient_3_MDSsAML")

for (file in scRNA_file){
        seurat_data <- Read10X(data.dir = paste0(scRNA_dir, file))
        seurat_obj <- CreateSeuratObject(counts = seurat_data$'Gene Expression', 
                                         min.features = 100, 
                                         project = file)
        assign(file, seurat_obj)
}

#merge all the Seuratobject 
merged_seurat <- merge(x = Patient_17_MDS, 
                       y = c(Patient_17_sAML, Patient_3_MDSsAML),
                       add.cell.id = c("Patient_17_MDS", "Patient_17_sAML", "Patient_3_MDSsAML"))
					   

##to see metadata created automatically by Seurant 
#View(test@meta.data)

####################################
##filter low quality cells###  name of seurat object == test 
####################################

########
#1. Add quality metrics 
########
##Add Novelty score to metadata slot 
#test[["log10GenesPerUMI"]] <- log10(test$nFeature_RNA) / log10(test$nCount_RNA)
#alternatve way to reach to metadata slot
test$log10GenesPerUMI <- log10(test$nFeature_RNA) / log10(test$nCount_RNA)

##Add Mitochondrial Gene percentage to metadata slot 
test$mitoRatio <- PercentageFeatureSet(object = test, pattern = "^MT-")
#test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^MT-")
#VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

########
#2. create metadata 
########

##separate metadata from seurat object
metadata <- test@meta.data
metadata$cells <- rownames(metadata)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^Patient_17_MDS"))] <- "Patient_17_MDS" ##this does not work when you did not merge samples 
#metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

## add metadata back to seurat object
test@meta.data <- metadata

#######
#3. visualise the quality metrics
######


######
#4. filtering (cells and genes)  
######
##filter out cells 
filtered_test <- subset(x = test, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI >= 0.80) & 
                           (mitoRatio <= 15))	   
#test <- subset(test, subset = nFeature_RNA > 250 & nCount_RNA > 500 & log10GenesPerUMI > 0.8 & percent.mt < 0.15)

##filter out genes
counts <- GetAssayData(object = filtered_test, slot = "counts")  #get count data from seurat object
nonzero <- counts > 0  #transform into logical matrix (values are . or | instead of number/quantity)
keep_genes <- Matrix::rowSums(nonzero) >= 10 #take genes with more than 10 cells which express the gene and returns TRUE FALSE matrix
filtered_counts <- counts[keep_genes, ]
filtered_test <- CreateSeuratObject(filtered_counts, meta.data = filtered_test@meta.data)

