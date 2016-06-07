#include <stdio.h>
#include <stdlib.h>

#include "../common.h"

static uint32_t ecc_compute_overall_parity_x(matrix_entry element)
{
  uint32_t *data = (uint32_t*)&element;
  return __builtin_parity(data[0] ^ data[1] ^ data[2] ^ data[3]);
}

#define TABLE_BITS 8
#define TABLE_SIZE (1<<TABLE_BITS)
static uint8_t table[TABLE_SIZE];
/*
{
  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,
  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,

  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,
  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,

  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,
  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,

  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,
  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0,
 };
*/


// Initialize ECC for a sparse matrix
void init_matrix_ecc(sparse_matrix M)
{
  for (int i = 0; i < TABLE_SIZE; i++)
  {
    table[i] = __builtin_parity(i);
  }

  // Add ECC protection to matrix elements
  for (unsigned i = 0; i < M.nnz; i++)
  {
    matrix_entry element = M.elements[i];

    // Compute overall parity bit for whole codeword
    element.col |= ecc_compute_overall_parity_x(element) << 31;

    M.elements[i] = element;
  }
}

// Sparse matrix vector product
// Multiplies `matrix` by `vector` and stores answer in `result`
// The matrix and vector dimensions are `N`
void spmv_x(sparse_matrix matrix, double *vector, double *result, unsigned N)
{
  // Initialize result vector to zero
  for (unsigned i = 0; i < N; i++)
    result[i] = 0.0;

  asm goto(
    "add     r5, %[ep], %[nnz], lsl #4\n"

    ".LOOP_BODY:\n\t"

    "ldr     r2, [%[ep]]\n\t"
    "ldr     r1, [%[ep], #4]\n\t"

    "eor     r0, r2, r1\n\t"
    "ldr     r4, [%[ep], #8]\n\t"
    "eor     r0, r0, r4\n\t"
    "ldr     r4, [%[ep], #12]\n\t"
    "eor     r0, r0, r4\n\t"

#if 1
    "eor     r0, r0, r0, lsr #16\n\t"
    "eor     r0, r0, r0, lsr #8\n\t"
    //"eor     r0, r0, r0, lsr #4\n\t"
    //"eor     r0, r0, r0, lsr #2\n\t"
    //"eor     r0, r0, r0, lsr #1\n\t"
    //"and     r0, r0, #0x1\n\t"
    "and     r0, r0, #0xFF\n\t"
    "ldrB    r0, [%[table], r0]\n\t"
#else
    "bl      __paritysi2\n\t"
#endif

    "cbnz    r0, %l[ERROR]\n\t"

    "and     r2, r2, #0x00FFFFFF\n\t"

    "add     r2, %[result], r2, lsl #3\n\t"
    "add     r1, %[vector], r1, lsl #3\n\t"
    "vldr.64 d5, [r2]\n\t"
    "vldr.64 d6, [%[ep], #8]\n\t"
    "vldr.64 d7, [r1]\n\t"
    "vmla.f64        d5, d6, d7\n\t"
    "vstr.64 d5, [r2]\n\t"

    "add     %[ep], %[ep], #16\n\t"
    "cmp     %[ep], r5\n\t"
    "bne     .LOOP_BODY\n"
    :
    : [result] "r" (result), [table] "r" (table),
      [ep] "r" (matrix.elements),
      [vector] "r" (vector),
      [nnz] "r" (matrix.nnz)
    : "cc", "r0", "r1", "r2", "r3", "r4", "r5", "lr", "d5", "d6", "d7", "memory"
    : ERROR
    );

  return;

ERROR:
  printf("[ECC] error detected at index (something)\n");
  exit(1);

  // Check overall parity bit
  //if (ecc_compute_overall_parity_x(element))
  //{
  //  printf("[ECC] error detected at index %d\n", i);
  //  exit(1);
  //}

  // Mask out ECC from high order column bits
  //element.col &= 0x00FFFFFF;
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
  int i = 0;
  do
  {
    // Load non-zero element
    matrix_entry element = matrix.elements[i];

    // Check overall parity bit
    //if (ecc_compute_overall_parity_x(element))

    uint32_t x = element.col ^ element.row;
    x ^= ((uint32_t*)(&element.value))[0];
    x ^= ((uint32_t*)(&element.value))[1];
    //if (__builtin_parity(x))
    x ^= x>>16;
    x ^= x>>8;
    //x ^= x>>4;
    //x ^= x>>2;
    //x ^= x>>1;
    //x &= 1;
    x &= 0xFF;
    if (table[x])
    {
      goto error;
    }

    // Mask out ECC from high order column bits
    element.col &= 0x00FFFFFF;

    // Multiply element value by the corresponding vector value
    // and accumulate into result vector
    result[element.col] += element.value * vector[element.row];
  }
  while (++i < matrix.nnz);

  return;

error:
  printf("[ECC] error detected at index (something)\n");
  exit(1);
}
