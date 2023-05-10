
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


# HDF5 - random

#Create HDF5 data
# loading "1M_neurons_data.h5ad" obtained from python script "scanpytenx.py"
setwd("/mnt/spca/run_spca/M1/mem")
input <- commandArgs(trailingOnly = TRUE)

here()
#### 100k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_100k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_100k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M1_Random", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)
time.end <- proc.time()
time100k_random<- time.end - time.start
time100k_random


# elapsed time in minute
time100k_random[3]/60
head(random_pca$x[,1:2])


#### 500k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_500k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_500k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M1_Random", out_name)), append = FALSE, memory.profiling = TRUE)


invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time500k_random <- time.end - time.start
time500k_random
# elapsed time in minute
time500k_random[3]/60
head(random_pca$x[,1:2])


#### 1M


time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1000k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1000k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M1_Random", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time1000k_random <- time.end - time.start
time1000k_random

# elapsed time in minute
time1000k_random[3]/60
head(random_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1.3M"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1.3M", "_", "ila",".out")

Rprof(filename = here("output",paste0("M1_Random", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time1.3M_random <- time.end - time.start
time1.3M_random

# elapsed time in minute
time1.3M_random[3]/60
head(random_pca$x[,1:2])


# HDF5 - exact
#### 100k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_100k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_100k", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Exact", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time100k_exact<- time.end - time.start
time100k_exact

# elapsed time in minute
time100k_exact[3]/60
head(random_pca$x[,1:2])


#### 500k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_500k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_500k", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Exact", out_name)), append = FALSE, memory.profiling = TRUE)
invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time500k_exact <- time.end - time.start
time500k_exact
# elapsed time in minute
time500k_exact[3]/60
head(random_pca$x[,1:2])


#### 1M


time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1000k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1000k", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Exact", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time1000k_exact <- time.end - time.start
time1000k_exact

# elapsed time in minute
time1000k_exact[3]/60
head(random_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1.3M"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1.3M", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Exact", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)
time.end <- proc.time()
time1.3M_exact <- time.end - time.start
time1.3M_exact

# elapsed time in minute
time1.3M_exact[3]/60
head(random_pca$x[,1:2])


#HDF5 - irlba

### 100k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_100k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_100k", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time100k_irlba<- time.end - time.start
time100k_irlba

# elapsed time in minute
time100k_irlba[3]/60
head(random_pca$x[,1:2])


#### 500k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_500k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_500k", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time500k_irlba <- time.end - time.start
time500k_irlba
# elapsed time in minute
time500k_irlba[3]/60
head(random_pca$x[,1:2])


#### 1M


time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1000k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1000k", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time1000k_irlba <- time.end - time.start
time1000k_irlba

# elapsed time in minute
time1000k_irlba[3]/60
head(random_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1.3M"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input,"_1.3M", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)
time.end <- proc.time()
time1.3M_irlba <- time.end - time.start
time1.3M_irlba

# elapsed time in minute
time1.3M_irlba[3]/60
head(random_pca$x[,1:2])



data100_random<- data.frame(dataset = "TENxBrain_100k",
                            ncells = 100000,
                            method = "random",
                            elapsed_time = time100k_random[3])

data500_random<- data.frame(dataset = "TENxBrain_500k",
                            ncells = 500000,
                            method = "random",
                            elapsed_time = time500k_random[3])

data1000_random<- data.frame(dataset = "TENxBrain_1000k",
                             ncells = 1000000,
                             method = "random",
                             elapsed_time = time1000k_random[3])

data1300_random<- data.frame(dataset = "TENxBrain_1.3M",
                             ncells = 1300000,
                             method = "random",
                             elapsed_time = time1.3M_random[3])

data_random <- rbind(data100_random, data500_random, data1000_random, data1300_random)
data_random

data100_exact<- data.frame(dataset = "TENxBrain_100k",
                           ncells = 100000,
                           method = "exact",
                           elapsed_time = time100k_exact[3])

data500_exact<- data.frame(dataset = "TENxBrain_500k",
                           ncells = 500000,
                           method = "exact",
                           elapsed_time = time500k_exact[3])

data1000_exact<- data.frame(dataset = "TENxBrain_1000k",
                            ncells = 1000000,
                            method = "exact",
                            elapsed_time = time1000k_exact[3])

data1300_exact<- data.frame(dataset = "TENxBrain_1.3M",
                            ncells = 1300000,
                            method = "exact",
                            elapsed_time = time1.3M_exact[3])

data_exact <- rbind(data100_exact, data500_exact, data1000_exact, data1300_exact)
data_exact 


data100_irlba<- data.frame(dataset = "TENxBrain_100k",
                           ncells = 100000,
                           method = "irlba",
                           elapsed_time = time100k_irlba[3])

data500_irlba<- data.frame(dataset = "TENxBrain_500k",
                           ncells = 500000,
                           method = "irlba",
                           elapsed_time = time500k_irlba[3])

data1000_irlba<- data.frame(dataset = "TENxBrain_1000k",
                            ncells = 1000000,
                            method = "irlba",
                            elapsed_time = time1000k_irlba[3])

data1300_irlba<- data.frame(dataset = "TENxBrain_1.3M",
                            ncells = 1300000,
                            method = "irlba",
                            elapsed_time = time1.3M_irlba[3])

data_irlba <- rbind(data100_irlba, data500_irlba, data1000_irlba, data1300_irlba)
data_irlba

data_M1 <- rbind(data_random, data_exact, data_irlba)
data_M1




#write.table(data_M1, file = paste("time_table", input, ".txt"))



