#include "../common.h"

// Initialize ECC for a sparse matrix
void init_matrix_ecc(sparse_matrix M)
{
  // Not using ECC - nothing to do
}

// Sparse matrix vector product
// Multiplies `matrix` by `vector` and stores answer in `result`
void spmv(sparse_matrix matrix, double *vector, double *result)
{
  // Initialize result vector to zero
  for (unsigned i = 0; i < matrix.N; i++)
    result[i] = 0.0;

  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < matrix.nnz; i++)
  {
    // Load non-zero element
    matrix_entry element = matrix.elements[i];

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result[element.col] += element.value * vector[element.row];
  }
}
