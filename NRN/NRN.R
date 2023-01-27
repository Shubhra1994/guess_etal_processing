##################################################
## Functions to performed nearest rank neighbors (NRN)
##################################################
#' 
#' @description NRN provided a appraoch to project  FACS produced data in to ABseq data.
#' For this the expression of each antibody was separately normalized both in the FACS and 
#' in the Abseq dataset using a rank-based approach. In particular, sample ranks were computed
#' and divided by the total number of samples, i.e. data was mapped to a scale from 0 to 1 
#' where 0 indicates lowest expression within the dataset, and 1 indicates highest expression. 
#' Within this normalized gene expression space, the cosine distance between any cell from the
#' Abseq (reference) dataset and the FC (query) dataset was computed, the four nearest reference
#' neighbours of every query cell were identified and the average position of these neighbours 
#' in uMAP and pseudotime space was computed using scmap.Subsequently, the average Euclidean 
#' distance of the reference neighbours in MOFA space was computed to identify cells with inconsistent
#'  mapping results. These cells can be later removed based on a user-defined threshold.
#' 
#' This function contains all steps required for performin NRN.
#' 
#' @details Please have a look at the README file for an example.
#' @name NRN
#' 
#' 
#' 
library(tidyr)
library(ggplot2)
library(SingleCellExperiment)
library(scmap)

#' @title projectOnUmap
#' @name projectOnUmap
#' @description Project FACS data into ABseq dataset 
#' @param FACS FACS data  as a dataframe
#' @param abseqSeurat a Seurat object of the ABseq dataset
#' @param pst  dataframe with pseudotime values per cell or any other variable to calculate
#' @param cd34_positivity_threshold  Threshold for CD34-AB
#' @param nnn  number of neirest neigbouhrs to take into account
#' @param k  number of clusters per group for k-means clusterin
#' @return Dataframe containing the calculated coordinates and pseudotime or variable feed as pst
#' @export
#' @examples 
##See read me
projectOnUmap <- function(FACS, abseqSeurat, pst=NULL, cd34_positivity_threshold = 1.7, nnn=4, k= 8500) {
  
  
  facsSeurat <- CreateSeuratObject(t(FACS), assay = "AB")
  
  quantile_normalisation <- function(df){
    df_rank <- apply(df,2,rank,ties.method="min")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)
    
    index_to_mean <- function(my_index, my_mean){
      return(my_mean[my_index])
    }
    
    df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
    rownames(df_final) <- rownames(df)
    return(df_final)
  }
  
  df_rank <- apply(t(facsSeurat@assays$AB@data),2,rank,ties.method="min")
  
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x))) }
  
  facsSeurat_AB <- t(apply(t(facsSeurat@assays$AB@data),2,rank,ties.method="min"))
  facsSeurat_AB <-t(apply(facsSeurat_AB,1, normalize))
  rownames(facsSeurat_AB)<-rownames(facsSeurat@assays$AB@data)
  colnames(facsSeurat_AB)<-colnames(facsSeurat@assays$AB@data)
  
  ##Subset Data to Progenitors
  toeval <- sprintf("abseqSeurat_CD34Pos<-subset(abseqSeurat,subset=`CD34-AB`>%f)",cd34_positivity_threshold)
  eval(str2expression(toeval))
  
  
  All_AB <- t(apply(t(abseqSeurat_CD34Pos@assays$AB@data[rownames(abseqSeurat_CD34Pos@assays$AB@data) %in% rownames(facsSeurat@assays$AB@data),]),2,rank,ties.method="min"))
  All_AB <-t(apply(All_AB,1, normalize))
  rownames(All_AB)<-rownames(abseqSeurat_CD34Pos@assays$AB@data[rownames(abseqSeurat_CD34Pos@assays$AB@data) %in% rownames(facsSeurat@assays$AB@data),])
  colnames(All_AB)<-colnames(abseqSeurat_CD34Pos@assays$AB@data[rownames(abseqSeurat_CD34Pos@assays$AB@data) %in% rownames(facsSeurat@assays$AB@data),])
  
  
  
  facsSeurat_AB_2<-as.matrix(facsSeurat_AB[rownames(All_AB),])
  rownames(facsSeurat_AB_2)<-paste(rownames(facsSeurat_AB_2),"_2")
  
  sce_Culture <- SingleCellExperiment(assays = list(normcounts = rbind(facsSeurat_AB_2,as.matrix(facsSeurat_AB[rownames(All_AB),]))))
  logcounts(sce_Culture) <- normcounts(sce_Culture)
  
  rowData(sce_Culture)$feature_symbol <- rownames(sce_Culture)
  
  sce_Culture <- sce_Culture[!duplicated(rownames(sce_Culture)), ]
  
  All_AB_2<-as.matrix(All_AB)
  rownames(All_AB_2)<-paste(rownames(All_AB_2),"_2")
  sce_All <- SingleCellExperiment(assays = list(normcounts = rbind(as.matrix(All_AB),All_AB_2)))
  logcounts(sce_All) <- normcounts(sce_All)
  # use gene names as feature symbols
  rowData(sce_All)$feature_symbol <- rownames(sce_All)
  # remove features with duplicated names
  sce_All <- sce_All[!duplicated(rownames(sce_All)), ]
  
  sce_Culture<-setFeatures(sce_Culture,features =  rownames(sce_Culture))
  sce_Culture <- indexCell(sce_Culture)
  
  sce_All<-setFeatures(sce_All,features =  rownames(sce_All))
  sce_All <- indexCell(sce_All,k=k)
  
  
  Culture_Map <- scmapCell(
    projection = sce_Culture,
    index_list = list(
      sce_All = metadata(sce_All)$scmap_cell_index
    ),
    w = nnn)
  
  #check for the nearest neighbors in FACS space that they are consistent; consistent meaning: small distances
  dd <- dist(Embeddings(abseqSeurat_CD34Pos, reduction = "MOFA"))
  
  dist.neighbors <- apply(Culture_Map$sce_All$cells,2,function(ii) {
    n <- attr(dd, "Size")
    subdist <- sapply(ii, function(i) {
      sapply(ii, function(j) {
        if (i < j) dd[n*(i-1) - i*(i-1)/2 + j-i] else NA
      })
    })
    mean(subdist, na.rm=T)
  })
  
  #calculate Coordinates
  Calc<-function(id,cult,vec){
    x=mean(vec[cult[,id]])
    return(x)
  }
  
  
  UmapX <- sapply(colnames(Culture_Map$sce_All[[1]]), function(id) {
    Calc(id,Culture_Map$sce_All[[1]], abseqSeurat_CD34Pos@reductions$MOFAUMAP@cell.embeddings[,1])
  })
  
  
  UmapY <- sapply(colnames(Culture_Map$sce_All[[1]]), function(id) {
    Calc(id,Culture_Map$sce_All[[1]], abseqSeurat_CD34Pos@reductions$MOFAUMAP@cell.embeddings[,2])
  })
  
  if (!is.null(pst)) {
    pst <- pst[colnames(abseqSeurat_CD34Pos),]
    pst.projected <- t(sapply(colnames(Culture_Map$sce_All[[1]]), function(id) {
      apply(pst, 2, function(xxx) Calc(id,Culture_Map$sce_All[[1]],xxx))
    }))
    CultureUMAP<-data.frame(UMAPX=UmapX,UMAPY=UmapY,cell=colnames(Culture_Map$sce_All[[1]]), celltype_C=paste0("c",Idents(facsSeurat)),orig="FACS",celltype=paste0("c",Idents(facsSeurat)), pst.projected,dist.neighbors = dist.neighbors)
    HealthyUMAP<-data.frame(UMAPX=abseqSeurat_CD34Pos@reductions$MOFAUMAP@cell.embeddings[,1],UMAPY=abseqSeurat_CD34Pos@reductions$MOFAUMAP@cell.embeddings[,2],
                            cell=rownames(abseqSeurat_CD34Pos@reductions$MOFAUMAP@cell.embeddings), celltype=Idents(abseqSeurat_CD34Pos),orig="Abseq",celltype_C=NA,pst,dist.neighbors = NA)
    
  } else {
    CultureUMAP<-data.frame(UMAPX=UmapX,UMAPY=UmapY,cell=colnames(Culture_Map$sce_All[[1]]), celltype_C=paste0("c",Idents(facsSeurat)),orig="FACS",celltype=paste0("c",Idents(facsSeurat)), dist.neighbors = dist.neighbors)
    HealthyUMAP<-data.frame(UMAPX=abseqSeurat_CD34Pos@reductions$MOFAUMAP@cell.embeddings[,1],UMAPY=abseqSeurat_CD34Pos@reductions$MOFAUMAP@cell.embeddings[,2],
                            cell=rownames(abseqSeurat_CD34Pos@reductions$MOFAUMAP@cell.embeddings), celltype=Idents(abseqSeurat_CD34Pos),orig="Abseq",celltype_C=NA, dist.neighbors = NA)
    
  }
  
  
  
  projected <- rbind(CultureUMAP, HealthyUMAP)
  return(projected)  
  
}
