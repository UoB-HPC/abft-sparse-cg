#include <stdio.h>
#include <stdlib.h>

#include "../common.h"
#include "../ecc.h"

// Initialize ECC for a sparse matrix
void init_matrix_ecc(sparse_matrix M)
{
  // Add ECC protection to matrix elements
  for (unsigned i = 0; i < M.nnz; i++)
  {
    csr_colval colval;
    colval.value = M.values[i];
    colval.column = M.cols[i];

    // Compute overall parity bit for whole codeword
    colval.column |= ecc_compute_overall_parity(colval) << 31;

    M.cols[i] = colval.column;
  }
}

// Sparse matrix vector product
// Multiplies `matrix` by `vector` and stores answer in `result`
// The matrix and vector dimensions are `N`
void spmv(sparse_matrix matrix, double *vector, double *result, unsigned N)
{
  for (unsigned row = 0; row < N; row++)
  {
    double tmp = 0.0;

    uint32_t start = matrix.rows[row];
    uint32_t end   = matrix.rows[row+1];
    for (uint32_t i = start; i < end; i++)
    {
      csr_colval colval;
      colval.value = matrix.values[i];
      colval.column = matrix.cols[i];

      // Check overall parity bit
      if (ecc_compute_overall_parity(colval))
      {
        printf("[ECC] error detected at index %d\n", i);
        exit(1);
      }

      // Mask out ECC from high order column bits
      colval.column &= 0x00FFFFFF;

      tmp += matrix.values[i] * vector[colval.column];
    }

    result[row] = tmp;
  }
}
