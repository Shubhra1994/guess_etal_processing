# NRN

NRN provides an appraoch to project FACS produced data in to ABseq data. For this the expression of each antibody was separately normalized both in the FACS and  in the Abseq dataset using a rank-based approach. In particular, sample ranks were computed and divided by the total number of samples, i.e. data was mapped to a scale from 0 to 1  where 0 indicates lowest expression within the dataset, and 1 indicates highest expression.  Within this normalized gene expression space, the cosine distance between any cell from the Abseq (reference) dataset and the FC (query) dataset was computed, the four nearest reference
neighbours of every query cell were identified and the average position of these neighbours 
in uMAP and pseudotime space was computed using scmap.Subsequently, the average Euclidean  distance of the reference neighbours in MOFA space was computed to identify cells with inconsistent  mapping results. These cells can be later removed based on a user-defined threshold.

## Quickstart
```R
Seurat_ABSEQ<-read.rds('path/to/your/ABseq_DATASET/folder')
facs <- read.csv('path/to/your/FACS_DATA/folder',stringsAsFactors = F)
rownames(facs) <- facs$CELL_NAME
projected <- projectOnUmap(norm.facs, Seurat_ABSEQ, pst = NULL, nnn=4, k =8000) #

```

