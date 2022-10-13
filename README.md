# MatrixTools

R function utilities to manipulate sparse matrices

### as_matrix
Convert a sparse matrix to a dense matrix. Use this when the as.matrix function throws an error of "Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105"

### dgT2dgCMatrix
Convert a dgTMatrix to a dgCMatrix.

### cbind_dgCMatrices
Merge (cbind) two dgCmatrix objects.
