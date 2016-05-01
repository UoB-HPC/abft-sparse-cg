#include <stdio.h>
#include <stdlib.h>

#include "common.h"

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
  for (unsigned i = 0; i < N; i++)
    result[i] = 0.0;

  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < matrix.nnz; i++)
  {
    // Load non-zero element
    matrix_entry element = matrix.elements[i];

    // Check index order constraints
    // Skip last row
    if (i < matrix.nnz - 1)
    {
      // Compare this row index to the next row index
      uint32_t next_row = matrix.elements[i+1].row;
      if (element.row > next_row)
      {
        printf("row index order violated for (%d,%d)\n",
               element.col, element.row);
        exit(1);
      }
      else if (element.row == next_row)
      {
        // Compare this column index to the next column index (if same row)
        uint32_t next_col = matrix.elements[i+1].col;
        if (element.col >= next_col)
        {
          printf("column index order violated for (%d,%d)\n",
                 element.col, element.row);
          exit(1);
        }
      }
    }

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result[element.col] += element.value * vector[element.row];
  }
}
