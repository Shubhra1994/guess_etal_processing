#every single cell experiment will generate: barcodes.tsv: cellular barcodes present in dataset
#features.tsv: IDs of quantified genes
#matrix.mtx: a matrix of count values, where rows are associated with the gene IDs above and columns correspond to the cellular barcodes. Note that there are many zero values in this matrix.

library("Matrix")
library("tidyverse")

# Read in `matrix.mtx`
counts <- readMM("GSM6229643_Patient_3_MDS_and_Patient_3_sAML_matrix.mtx.gz")
# Read in `genes.tsv`
genes <- read_tsv("GSM6229643_Patient_3_MDS_and_Patient_3_sAML_features.tsv.gz", col_names = FALSE)
gene_ids <- genes$X1
# Read in `barcodes.tsv`
cell_ids <- read_tsv("GSM6229643_Patient_3_MDS_and_Patient_3_sAML_barcodes.tsv.gz", col_names = FALSE)$X1


# Make the column names as the cell IDs and the row names as the gene IDs
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids

#alternate way to read the mtx file using seurat internal; the name reading is hardcoded in 
data_dir <- '/mnt/beegfs/sbhattacharya/scRNA_standardization/Guess_copy_etal/Patient_3'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir,gene.column = 2,  cell.column = 1,  strip.suffix = FALSE)
seurat_object = CreateSeuratObject(counts = data$`Gene Expression`) #data[[1]],data[[2]]
seurat_object[['RNA']] = CreateAssayObject(counts = data$`Gene Expression`)
#seurat_object2[['Ab']] = CreateAssayObject(counts = data$`Antibody Capture`)
#seurat_object2 = CreateSeuratObject(counts = data$`Antibody Capture`)
