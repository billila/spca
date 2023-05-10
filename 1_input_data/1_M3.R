# 1_M3

setwd("~/")

library(TENxBrainData)
tenx <- TENxBrainData()

rownames(tenx)<- rowData(tenx)$Symbol

# load the data preprocessed in python to filter the matrix
dati <- "1M_neurons_data.h5ad"
dati <- readH5AD(dati, use_hdf5 = TRUE)

tenx <- tenx[rownames(dati),]
assay(tenx)

# normalization 
tenx <- logNormCounts(tenx)

y <- counts(tenx)

# save y as TileDB
library(TileDBArray)

tdb.mat <- as(y, "TileDBArray")

# save as sparse counts matrix
writeTileDBArray(y, path= "tileDB_1.3M.tdb", sparse = TRUE)


ncells <- c(100, 500, 1000)
cellidx <- lapply(ncells, function(n) sample(colnames(y), n * 1000))
sapply(cellidx, length)

sub100k <- y[, cellidx[[1]]]
writeTileDBArray(sub100k, path= "tileDB_100k.tdb", sparse = TRUE)

sub500k <- y[, cellidx[[2]]]
writeTileDBArray(sub500k, path= "tileDB_500k.tdb", sparse = TRUE)

sub1M <- y[, cellidx[[3]]]
writeTileDBArray(sub1M, path= "tileDB_1M.tdb", sparse = TRUE)

# example: load 100k matrix
tenxDB <- TileDBArray("~/tileDB_100k.tdb")

# check the class
class(tenxDB)
# check the sparsity
is_sparse(tenxDB)

# load the TileDB object to compute PCA

setwd("/mnt/spca/M3/time")
mat100k <- TileDBArray("tileDB_100k.tdb")
mat500k <- TileDBArray("tileDB_500k.tdb")
mat1M <- TileDBArray("tileDB_1M.tdb")
mat <- TileDBArray("tileDB_1.3M.tdb")