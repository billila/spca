# time M4
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

# load TENxBrainData
tenx <- TENxMatrix("/mnt/spca/M4/time/1M_neurons_filtered_gene_bc_matrices_h5.h5", group="mm10")


# transform tenx in a sce S4 object
library(SingleCellExperiment)

sce <- SingleCellExperiment(assays = list (counts = tenx),
                            colData = tenx@seed@dimnames[[2]],
                            rowData = tenx@seed@dimnames[[1]])

setwd("~/")
dati <- "1M_neurons_data.h5ad"
dati <- readH5AD(dati, use_hdf5 = TRUE)

rownames(dati)<- rowData(dati)$gene_ids
sce <- sce[rownames(dati),]

sce <- logNormCounts(sce)
y <- counts(sce)

# save y as TileDB
library(TileDBArray)

tdb.mat <- as(y, "TileDBArray")

# save as sparse counts matrix
#setwd("/mnt/spca/M4/time")
setwd("/home/ilaria/Documents/Scalable PCA/M4/time")
writeTileDBArray(y, path= "/home/ilaria/Documents/Scalable PCA/M4/time/tileDB_1.3M.tdb", sparse = TRUE)


ncells <- c(100, 500, 1000)
cellidx <- lapply(ncells, function(n) sample(colnames(y), n * 1000))
sapply(cellidx, length)

sub100k <- y[, cellidx[[1]]]
writeTileDBArray(sub100k, path= "/home/ilaria/Documents/Scalable PCA/M4/time/tileDB_100k.tdb", sparse = TRUE)

sub500k <- y[, cellidx[[2]]]
writeTileDBArray(sub500k, path= "/home/ilaria/Documents/Scalable PCA/M4/time/tileDB_500k.tdb", sparse = TRUE)

sub1M <- y[, cellidx[[3]]]
writeTileDBArray(sub1M, path= "/home/ilaria/Documents/Scalable PCA/M4/time/tileDB_1M.tdb", sparse = TRUE)


#set your directory
tenxDB <- TileDBArray("~/tileDB_100k.tdb")

# check the class
class(tenxDB)
# check the sparsity
is_sparse(tenxDB)



