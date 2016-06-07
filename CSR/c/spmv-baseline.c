#include "../common.h"

// Initialize ECC for a sparse matrix
void init_matrix_ecc(sparse_matrix M)
{
  // Not using ECC - nothing to do
}

// Sparse matrix vector product
// Multiplies `matrix` by `vector` and stores answer in `result`
// The matrix and vector dimensions are `N`
void spmv(sparse_matrix matrix, double *vector, double *result, unsigned N)
{
  // Initialize result vector to zero
  //for (unsigned i = 0; i < N; i++)
  //  result[i] = 0.0;

  for (unsigned row = 0; row < N; row++)
  {
    double tmp = 0.0;

    uint32_t start = matrix.rows[row];
    uint32_t end   = matrix.rows[row+1];
    for (uint32_t i = start; i < end; i++)
    {
      uint32_t col = matrix.cols[i];
      tmp += matrix.values[i] * vector[col];
    }

    result[row] = tmp;
  }

  // // Loop over non-zeros in matrix
  // for (unsigned i = 0; i < matrix.nnz; i++)
  // {
  //   // Load non-zero element
  //   matrix_entry element = matrix.elements[i];
  //
  //   // Multiply element value by the corresponding vector value
  //   // and accumulate into result vector
  //   result[element.col] += element.value * vector[element.row];
  // }
}
