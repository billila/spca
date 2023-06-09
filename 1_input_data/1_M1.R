# 1_M1
library(zellkonverter)
# input M1

setwd("~/")

library(TENxBrainData)
tenx <- TENxBrainData()

rownames(tenx)<- rowData(tenx)$Symbol

# load the data preprocessed in python to filter the matrix
dati <- "1M_neurons_data.h5ad"
dati <- readH5AD(dati, use_hdf5 = TRUE)

tenx <- tenx[rownames(dati),]
SummarizedExperiment::assay(tenx)

# normalization 
tenx <- logNormCounts(tenx)


if(!file.exists(here("metodo2/data"))){
  dir.create(here("metodo2/data"))
}

## Check that the counts object is a HDF5Array
seed(logcounts(tenx))

# create subset
set.seed(138)
ncells <- c(100, 500, 1000, 1306.127) 
cellidx <- lapply(ncells, function(n) sample(colnames(tenx), n * 1000))
sapply(cellidx, length)



## save the subset as HDF5

if(!file.exists(here("metodo2/data/subset/TENxBrainDataSE"))) {
  dir.create(here("metodo2/datasubset/TENxBrainDataSE"), recursive = TRUE)
}

for(i in seq_along(cellidx)) {
  if(!file.exists(here("metodo2/data/subset/TENxBrainDataSE", paste0("TENxBrainData_", ncells[[i]], "k")))) {
    dir.create(here("metodo2/data/subset/TENxBrainDataSE", paste0("TENxBrainData_", ncells[[i]], "k")), recursive = TRUE)
  }
  
  sub <- tenx[, cellidx[[i]]]
  
  saveHDF5SummarizedExperiment(sub, 
                               dir = here("metodo2/data/subset/TENxBrainDataSE",
                                          paste0("TENxBrainData_", ncells[[i]], "k")), 
                               prefix="", replace=TRUE, 
                               #chunkdim=c(dim(counts(sub))[1],1), 
                               level=NULL, verbose=FALSE)
  
  
  
  print(ncells[[i]])
}


# load an Large singleCellExperiment object saved before with saveHDF5SummarizedExperiment

# 100k
sce_100k <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_100k"), prefix="")

# 500k
sce_500k <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_500k"), prefix="")

# 1M
sce_1M <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1000k"), prefix="")

# 1.3M
sce_1.3M <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1306.127k"), prefix="")
