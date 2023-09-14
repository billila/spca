import sklearn as sk
import h5py as h5
from sklearn.decomposition import IncrementalPCA
import time
import os
import sys
import scanpy as sc
from scipy import sparse

n_components = 50

adata = sc.read_csv("/mnt/spca/run_spca/data_csv/mat100k.csv", delimiter=',', first_column_names= True, dtype='float32')

X_sparse = sparse.csr_matrix(adata.X)

start_time = time.time()

op = sk.decomposition.PCA(n_components = n_components, svd_solver = "full")
X_transformed = op.fit_transform(X_sparse)
X_transformed.shape

sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


adata = sc.read_csv("/mnt/spca/run_spca/data_csv/mat500k.csv", delimiter=',', first_column_names= True, dtype='float32')

X_sparse = sparse.csr_matrix(adata.X)

start_time = time.time()

op = sk.decomposition.PCA(n_components = n_components, svd_solver = "full")
X_transformed = op.fit_transform(X_sparse)
X_transformed.shape

sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


adata = sc.read_csv("/mnt/spca/run_spca/data_csv/mat1M.csv", delimiter=',', first_column_names= True, dtype='float32')

X_sparse = sparse.csr_matrix(adata.X)

start_time = time.time()

op = sk.decomposition.PCA(n_components = n_components, svd_solver = "full")
X_transformed = op.fit_transform(X_sparse)
X_transformed.shape

sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


adata = sc.read_csv("/mnt/spca/run_spca/data_csv/mat13M.csv", delimiter=',', first_column_names= True, dtype='float32')

X_sparse = sparse.csr_matrix(adata.X)

start_time = time.time()

op = sk.decomposition.PCA(n_components = n_components, svd_solver = "full")
X_transformed = op.fit_transform(X_sparse)
X_transformed.shape

sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
