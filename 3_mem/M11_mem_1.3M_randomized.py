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


os.chdir("/mnt/spca/run_spca/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1.3M/")


start_time = time.time()

file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)
sa = sparse.csr_matrix(adata.X) 

# convert sparse matrix into anndata 
adata2 = ad.AnnData(sa)
sc.tl.pca(adata2, svd_solver = "randomized")
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
