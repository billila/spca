import sklearn as sk
import h5py as h5
from sklearn.decomposition import IncrementalPCA
from memory_profiler import memory_usage
from time import sleep

def f():
  dsname = "assay001"
  fn = "assays.h5"
  n_components = 50
  matref = h5.File(fn, mode = "r")
  #chunk_size = 10
  op = sk.decomposition.IncrementalPCA(n_components = n_components)
  #, batch_size = chunk_size
  op.fit_transform(matref[dsname])
  

mem_usage = memory_usage(f)
print('Memory usage (in chunks of .1 seconds): %s' % mem_usage)
print('Maximum memory usage: %s' % max(mem_usage))

