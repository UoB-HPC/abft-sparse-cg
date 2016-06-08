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
  for (unsigned row = 0; row < matrix.N; row++)
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
}
