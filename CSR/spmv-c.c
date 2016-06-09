#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "ecc.h"

static void spmv_baseline(sparse_matrix matrix, double *vector, double *result)
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

static void spmv_constraints(sparse_matrix matrix, double *vector, double *result)
{
  for (unsigned row = 0; row < matrix.N; row++)
  {
    double tmp = 0.0;

    uint32_t start = matrix.rows[row];
    uint32_t end   = matrix.rows[row+1];

    if (end > matrix.nnz)
    {
      printf("row size constraint violated for row %d\n", row);
      exit(1);
    }
    if (end < start)
    {
      printf("row order constraint violated for row%d\n", row);
      exit(1);
    }

    for (uint32_t i = start; i < end; i++)
    {
      uint32_t col = matrix.cols[i];

      if (col >= matrix.N)
      {
        printf("column size constraint violated at index %d\n", i);
        exit(1);
      }
      if (i < end-1)
      {
        if (matrix.cols[i+1] <= col)
        {
          printf("column order constraint violated at index %d\n", i);
          exit(1);
        }
      }

      tmp += matrix.values[i] * vector[col];
    }

    result[row] = tmp;
  }
}

static void spmv_sed(sparse_matrix matrix, double *vector, double *result)
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
        printf("[ECC] error detected at index %d\n", i);
        exit(1);
      }

      // Mask out ECC from high order column bits
      colval.column &= 0x00FFFFFF;

      tmp += colval.value * vector[colval.column];
    }

    result[row] = tmp;
  }
}

static void spmv_sec7(sparse_matrix matrix, double *vector, double *result)
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

      tmp += colval.value * vector[colval.column];
    }

    result[row] = tmp;
  }
}

static void spmv_sec8(sparse_matrix matrix, double *vector, double *result)
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

static void spmv_secded(sparse_matrix matrix, double *vector, double *result)
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

      // Check parity bits
      uint32_t overall_parity = ecc_compute_overall_parity(colval);
      uint32_t syndrome = ecc_compute_col8(colval);
      if (overall_parity)
      {
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
      colval.column &= 0x00FFFFFF;

      tmp += colval.value * vector[colval.column];
    }

    result[row] = tmp;
  }
}

void spmv(sparse_matrix matrix, double *vector, double *result)
{
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
