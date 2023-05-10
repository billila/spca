# 1_M8

# import library
import sklearn as sk
import h5py as h5
from sklearn.decomposition import IncrementalPCA

# you can find the "assays.h5" file in the file created
# with this function in M1: "saveHDF5SummarizedExperiment".

dsname = "assay001"
fn = "assays.h5"
matref = h5.File(fn, mode = "r")
