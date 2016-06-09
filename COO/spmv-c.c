#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "ecc.h"

static void spmv_baseline(sparse_matrix matrix, double *vector, double *result)
{
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

static void spmv_constraints(sparse_matrix matrix, double *vector, double *result)
{
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
        printf("row index order violated at index %d\n", i);
        exit(1);
      }
      else if (element.row == next_row)
      {
        // Compare this column index to the next column index (if same row)
        uint32_t next_col = matrix.elements[i+1].col;
        if (element.col >= next_col)
        {
          printf("column index order violated at index %d\n", i);
          exit(1);
        }
      }
    }

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result[element.col] += element.value * vector[element.row];
  }
}

static void spmv_sed(sparse_matrix matrix, double *vector, double *result)
{
  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < matrix.nnz; i++)
  {
    // Load non-zero element
    matrix_entry element = matrix.elements[i];

    // Check overall parity bit
    if (ecc_compute_overall_parity(element))
    {
      printf("[ECC] error detected at index %d\n", i);
      exit(1);
    }

    // Mask out ECC from high order column bits
    element.col &= 0x00FFFFFF;

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result[element.col] += element.value * vector[element.row];
  }
}

static void spmv_sec7(sparse_matrix matrix, double *vector, double *result)
{
  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < matrix.nnz; i++)
  {
    // Load non-zero element
    matrix_entry element = matrix.elements[i];

    // Check ECC
    uint32_t syndrome = ecc_compute_col8(element);
    if (syndrome)
    {
      // Unflip bit
      uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
      flip_bit(&element, bit);
      matrix.elements[i] = element;

      printf("[ECC] corrected bit %u at index %d\n", bit, i);
    }

    // Mask out ECC from high order column bits
    element.col &= 0x00FFFFFF;

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result[element.col] += element.value * vector[element.row];
  }
}

static void spmv_sec8(sparse_matrix matrix, double *vector, double *result)
{
  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < matrix.nnz; i++)
  {
    // Load non-zero element
    matrix_entry element = matrix.elements[i];

    // Check overall parity bit
    if (ecc_compute_overall_parity(element))
    {
      // Compute error syndrome from hamming bits
      uint32_t syndrome = ecc_compute_col8(element);
      if (syndrome)
      {
        // Unflip bit
        uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
        flip_bit(&element, bit);

        printf("[ECC] corrected bit %u at index %d\n", bit, i);
      }
      else
      {
        // Correct overall parity bit
        element.col ^= 0x1 << 24;

        printf("[ECC] corrected overall parity bit at index %d\n", i);
      }
      matrix.elements[i] = element;
    }

    // Mask out ECC from high order column bits
    element.col &= 0x00FFFFFF;

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result[element.col] += element.value * vector[element.row];
  }
}

static void spmv_secded(sparse_matrix matrix, double *vector, double *result)
{
  // Loop over non-zeros in matrix
  for (unsigned i = 0; i < matrix.nnz; i++)
  {
    // Load non-zero element
    matrix_entry element = matrix.elements[i];

    // Check parity bits
    uint32_t overall_parity = ecc_compute_overall_parity(element);
    uint32_t syndrome = ecc_compute_col8(element);
    if (overall_parity)
    {
      if (syndrome)
      {
        // Unflip bit
        uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
        flip_bit(&element, bit);

        printf("[ECC] corrected bit %u at index %d\n", bit, i);
      }
      else
      {
        // Correct overall parity bit
        element.col ^= 0x1 << 24;

        printf("[ECC] corrected overall parity bit at index %d\n", i);
      }
      matrix.elements[i] = element;
    }
    else
    {
      if (syndrome)
      {
        // Overall parity fine but error in syndrom
        // Must be double-bit error - cannot correct this
        printf("[ECC] double-bit error detected\n");
        exit(1);
      }
    }

    // Mask out ECC from high order column bits
    element.col &= 0x00FFFFFF;

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result[element.col] += element.value * vector[element.row];
  }
}

// Sparse matrix vector product
// Multiplies `matrix` by `vector` and stores answer in `result`
void spmv(sparse_matrix matrix, double *vector, double *result)
{
  // Initialize result vector to zero
  for (unsigned i = 0; i < matrix.N; i++)
    result[i] = 0.0;

  switch (matrix.mode)
  {
  case NONE:
    spmv_baseline(matrix, vector, result);
    break;
  case CONSTRAINTS:
    spmv_constraints(matrix, vector, result);
    break;
  case SED:
    spmv_sed(matrix, vector, result);
    break;
  case SEC7:
    spmv_sec7(matrix, vector, result);
    break;
  case SEC8:
    spmv_sec8(matrix, vector, result);
    break;
  case SECDED:
    spmv_secded(matrix, vector, result);
    break;
  }
}
