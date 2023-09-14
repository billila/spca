# M5_time


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


load("/mnt/spca/run_spca/M5_deferred/time/data_M5.RData")

#### 100k

time.start <- proc.time()


invisible(random_pca <- BiocSingular::runPCA(mat100k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time100k_random<- time.end - time.start
time100k_random

# elapsed time in minute
time100k_random[3]/60
head(random_pca$x[,1:2])


#### 500k

time.start <- proc.time()


invisible(random_pca <- BiocSingular::runPCA(mat500k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time500k_random <- time.end - time.start
time500k_random
# elapsed time in minute
time500k_random[3]/60
head(random_pca$x[,1:2])


#### 1M


time.start <- proc.time()


invisible(random_pca <- BiocSingular::runPCA(mat1M, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time1000k_random <- time.end - time.start
time1000k_random

# elapsed time in minute
time1000k_random[3]/60
head(random_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()


invisible(random_pca <- BiocSingular::runPCA(mat, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))
time.end <- proc.time()
time1.3M_random <- time.end - time.start
time1.3M_random

# elapsed time in minute
time1.3M_random[3]/60
head(random_pca$x[,1:2])


# HDF5 - exact
#### 100k

time.start <- proc.time()

invisible(random_pca <- BiocSingular::runPCA(mat100k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time100k_exact<- time.end - time.start
time100k_exact

# elapsed time in minute
time100k_exact[3]/60
head(random_pca$x[,1:2])


#### 500k

time.start <- proc.time()

invisible(random_pca <- BiocSingular::runPCA(mat500k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time500k_exact <- time.end - time.start
time500k_exact
# elapsed time in minute
time500k_exact[3]/60
head(random_pca$x[,1:2])


#### 1M


time.start <- proc.time()

invisible(random_pca <- BiocSingular::runPCA(mat1M, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time1000k_exact <- time.end - time.start
time1000k_exact

# elapsed time in minute
time1000k_exact[3]/60
head(random_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()

invisible(random_pca <- BiocSingular::runPCA(mat, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))
time.end <- proc.time()
time1.3M_exact <- time.end - time.start
time1.3M_exact

# elapsed time in minute
time1.3M_exact[3]/60
head(random_pca$x[,1:2])


#HDF5 - irlba

### 100k

time.start <- proc.time()


invisible(random_pca <- BiocSingular::runPCA(mat100k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time100k_irlba<- time.end - time.start
time100k_irlba

# elapsed time in minute
time100k_irlba[3]/60
head(random_pca$x[,1:2])


#### 500k

time.start <- proc.time()

invisible(random_pca <- BiocSingular::runPCA(mat500k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time500k_irlba <- time.end - time.start
time500k_irlba
# elapsed time in minute
time500k_irlba[3]/60
head(random_pca$x[,1:2])


#### 1M


time.start <- proc.time()

invisible(random_pca <- BiocSingular::runPCA(mat1M, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time1000k_irlba <- time.end - time.start
time1000k_irlba

# elapsed time in minute
time1000k_irlba[3]/60
head(random_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()

invisible(random_pca <- BiocSingular::runPCA(mat, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))
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

data_M5 <- rbind(data_random, data_exact, data_irlba)
data_M5



tmp <- commandArgs(trailingOnly = TRUE)
write.table(data_M5, file = tmp)



