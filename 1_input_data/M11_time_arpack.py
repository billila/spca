# M11 arpack

import numpy as np
import pandas as pd
import scanpy as sc
import sys
import time
import os
import numpy as np
from scipy import sparse
import anndata as ad

sc.settings.verbosity  = 3
sc.settings.set_figure_params(dpi = 70)
sc.logging.print_versions()


os.chdir("/mnt/spca/run_spca/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_100k/")






file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)
sa = sparse.csr_matrix(adata.X) 
adata2 = ad.AnnData(sa)

start_time = time.time()
sc.tl.pca(adata2, svd_solver = "arpack")
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


os.chdir("/mnt/spca/run_spca/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_500k/")



file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)
sa = sparse.csr_matrix(adata.X) 
adata2 = ad.AnnData(sa)

start_time = time.time()
sc.tl.pca(adata2, svd_solver = "arpack")
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


os.chdir("/mnt/spca/run_spca/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1000k/")


file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)
sa = sparse.csr_matrix(adata.X) 
adata2 = ad.AnnData(sa)

start_time = time.time()
sc.tl.pca(adata2, svd_solver = "arpack")
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


os.chdir("/mnt/spca/run_spca/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1.3M/")


file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)
sa = sparse.csr_matrix(adata.X) 
adata2 = ad.AnnData(sa)

start_time = time.time()
sc.tl.pca(adata2, svd_solver = "arpack")
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
