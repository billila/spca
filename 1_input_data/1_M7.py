# 1_M7

# import library
import numpy as np
import pandas as pd
import scanpy as sc

# setting graphics info
sc.settings.verbosity  = 3
sc.settings.set_figure_params(dpi = 70)
sc.logging.print_versions()

# you can find the "assays.h5" file in the file created
# with this function in M1: "saveHDF5SummarizedExperiment".
file = "assays.h5"
str = "assay001"

# read the file
adata=sc.read_hdf(file, str)
