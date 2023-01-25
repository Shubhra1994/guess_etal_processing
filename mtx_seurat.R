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
