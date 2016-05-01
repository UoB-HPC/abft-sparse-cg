#include "common.h"

// Initialize ECC for a sparse matrix
void init_matrix_ecc(sparse_matrix M)
{
  // Add ECC protection to matrix elements
  for (unsigned i = 0; i < M.nnz; i++)
  {
    matrix_entry element = M.elements[i];

    // TODO: Generate parity bits and store in high order column bits
    // Use the ecc_compute_col8 function to generate the parity bits
    // and then add them to the column index with the | operator (bitwise OR)

    M.elements[i] = element;
  }
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

    // TODO: Check ECC here
    // Use the ecc_compute_col8 function to check all of the parity bits
    // If there was an error use the ecc_correct_col8 function to
    // correct it

    // TODO: Mask out the parity bits from the high order column bits so that
    // we can use the column index
    // Use the & operator (bitwise AND) with an approriate mask

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result[element.col] += element.value * vector[element.row];
  }
}
