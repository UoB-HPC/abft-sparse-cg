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

    // Compute overall parity bit for whole codeword
    colval.column |= ecc_compute_overall_parity(colval) << 24;

    M.cols[i] = colval.column;
  }
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
      csr_colval colval;
      colval.value = matrix.values[i];
      colval.column = matrix.cols[i];

      // Check overall parity bit
      if (ecc_compute_overall_parity(colval))
      {
        // Compute error syndrome from hamming bits
        uint32_t syndrome = ecc_compute_col8(colval);
        if (syndrome)
        {
          // Unflip bit
          uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
          flip_bit(&colval, bit);

          printf("[ECC] corrected bit %u at index %d\n", bit, i);
        }
        else
        {
          // Correct overall parity bit
          colval.column ^= 0x1 << 24;

          printf("[ECC] corrected overall parity bit at index %d\n", i);
        }
        matrix.cols[i] = colval.column;
        matrix.values[i] = colval.value;
      }

      // Mask out ECC from high order column bits
      colval.column &= 0x00FFFFFF;

      tmp += colval.value * vector[colval.column];
    }

    result[row] = tmp;
  }
}
