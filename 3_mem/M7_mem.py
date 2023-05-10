import numpy as np
import pandas as pd
import scanpy as sc
from memory_profiler import memory_usage
from time import sleep


sc.settings.verbosity  = 3
sc.settings.set_figure_params(dpi = 70)
sc.logging.print_versions()


def f():
  file = "assays.h5"
  str = "assay001"
  sleep(.1)
  adata=sc.read_hdf(file, str)
  sleep(.1)
  sc.tl.pca(adata)
  sleep(.1)

mem_usage = memory_usage(f)
print('Memory usage (in chunks of .1 seconds): %s ' % mem_usage)
print('Maximum memory usage: %s ' % max(mem_usage))

with open('sample.txt', 'w') as f:
  f.write("%s = %s\n" %("dict", d))
  x = input()
