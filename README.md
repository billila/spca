# spca

With the advances in sequencing technology, the size and complexity of single-cell RNA-seq data are increasing, to the point that standard workflows are becoming too computationally demanding. In fact, datasets with millions of cells are becoming routine and they require workflows that operate out-of-memory.
However, existing software tools for single cells do not scale well to such large datasets. As an exemplary step, we consider dimensionality reduction, which is often one of the first steps of a typical analysis. Popular algorithms, such as principal component analysis (PCA), are based on singular value decomposition (SVD).
The classic implementations of SVD can be slow with large data matrices. Furthermore, they require the data to be loaded entirely into memory and therefore can be impossible to run with large datasets.

We tested the above methods on a real single-cell RNA-seq dataset from 10X Genomics that contains approximately 1.3 million cells and 30,000 genes isolated from the mouse brain.
We also created downsampled datasets (of sizes 100k, 500k, and 1M) to check the scalability of the methods (combination of input and algorithm) both in terms of time and memory consumption.

* The fold "0_data_preprocessing" contains the code to obtain the data used in this benchmark and the preprocessing.

* The fold "1_input_data" contains the code to obtain the input data for the PCA.

* The fold "2_time" contains the code to compute pca for each data type described in the fold "1_imput_data" and for each subsample of the data.

* The fold "3_time" contains the code to compute pca and the code to save the computational cost (GB).

  
![](https://github.com/billila/spca/blob/main/0_data_preprocessing/SCHEMA%20POSTER%20(1024%20%C3%97%20540%20px)%20(1024%20%C3%97%20640%20px)(4).png)

### LEGEND:
* M1 = bioc_hdf5
* M2 = bioc_dense
* M3 = TBD
* M4 = TBD
* M5 = bioc_sparse
* M5_deferred = bioc_sparse_deferred 
* M6 = bioc_sparse_array
* M6_deferred = bioc_sparse_array_deferred
* M7 = scanpy_dense
* M8 = sciktlearn_dense_IPCA
* M9 = TBD
* M10 = rapids_dense
* M11 = scanpy_sparse
* M12 = sciktlearn_sparse_IPCA
* M13 = rapids_sparse
* M14 = rspectra_dense
* M15 = rspectra_sparse
* M16 = sciktlearn_dense_exact
* M17 = sciktlearn_sparse_exact
