#include <stdio.h>

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

    // Generate ECC and store in high order column bits
    colval.column |= ecc_compute_col8(colval);

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

      // Check ECC
      uint32_t syndrome = ecc_compute_col8(colval);
      if (syndrome)
      {
        // Unflip bit
        uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
        flip_bit(&colval, bit);
        matrix.cols[i] = colval.column;
        matrix.values[i] = colval.value;

        printf("[ECC] corrected bit %u at index %d\n", bit, i);
      }

      // Mask out ECC from high order column bits
      colval.column &= 0x00FFFFFF;

      tmp += matrix.values[i] * vector[colval.column];
    }

    result[row] = tmp;
  }
}
