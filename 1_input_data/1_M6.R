# M6 PCA 

library(SingleCellExperiment)
library(zellkonverter)
library(BiocSingular)
library(here)
library(HDF5Array)
library(mbkmeans)
library(ggplot2)
library(scater)
library(scran)
library(BiocParallel)
library(DelayedMatrixStats)
library(rhdf5)
library(TileDBArray)
library(S4Arrays)
library(SparseArray)

# tileDB + BiocSingular 

# load TENxBrainData
tenx <- TENxMatrix("/mnt/spca/M4/time/1M_neurons_filtered_gene_bc_matrices_h5.h5", group="mm10")


# transform tenx in a sce S4 object
library(SingleCellExperiment)

sce <- SingleCellExperiment(assays = list (counts = tenx),
                            colData = tenx@seed@dimnames[[2]],
                            rowData = tenx@seed@dimnames[[1]])

setwd("/mnt/spca/M1/time/")
dati <- "1M_neurons_data.h5ad"
dati <- readH5AD(dati, use_hdf5 = TRUE)


rownames(dati)<- rowData(dati)$gene_ids
sce <- sce[rownames(dati),]

sce <- logNormCounts(sce)
y <- counts(sce)

mat <- as(y, "sparseMatrix")

mat <- SparseArray(mat)

ncells <- c(100, 500, 1000)
cellidx <- lapply(ncells, function(n) sample(colnames(mat), n * 1000))
sapply(cellidx, length)

mat100k <- mat[,cellidx[[1]]]
mat500k <- mat[,cellidx[[2]]]
mat1M <- mat[,cellidx[[3]]]

save.image("/mnt/spca/M6/time/data_M6.RData")
