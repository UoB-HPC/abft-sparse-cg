#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "ecc.h"

#define USE_PARITYSI2 0
#define USE_PARITY_TABLE 1

#if USE_PARITYSI2 + USE_PARITY_TABLE > 1
#error "Multiple implementations selected"
#endif

#if USE_PARITY_TABLE
static uint8_t PARITY_TABLE[256] =
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
#endif

static void spmv_baseline(sparse_matrix matrix, double *vector, double *result)
{
  printf("spmv_baseline not implemented\n");
  exit(1);
}

static void spmv_constraints(sparse_matrix matrix, double *vector, double *result)
{
  printf("spmv_constraints not implemented\n");
  exit(1);
}

static void __attribute__((noinline)) spmv_sed(sparse_matrix matrix, double *vector, double *result)
{
  asm goto(
    // Compute pointer to end of data
    "add     r5, %[ep], %[nnz], lsl #4\n"

    ".LOOP_BODY:\n\t"

    // Load column to r2 and row to r1
    "ldr     r2, [%[ep]]\n\t"
    "ldr     r1, [%[ep], #4]\n\t"

    // *** Parity check starts ***
    // Reduce data to 32-bits in r0
    "eor     r0, r2, r1\n\t"
    "ldr     r4, [%[ep], #8]\n\t"
    "eor     r0, r0, r4\n\t"
    "ldr     r4, [%[ep], #12]\n\t"
    "eor     r0, r0, r4\n\t"

#if USE_PARITYSI2
    // Use builtin to compute final parity
    "bl      __paritysi2\n\t"
#else
    "eor     r0, r0, r0, lsr #16\n\t"
    "eor     r0, r0, r0, lsr #8\n\t"

#if USE_PARITY_TABLE
    // Lookup final parity from table
    "and     r0, r0, #0xFF\n\t"
    "ldrB    r0, [%[PARITY_TABLE], r0]\n\t"
#else
    // Manually compute final parity
    "eor     r0, r0, r0, lsr #4\n\t"
    "eor     r0, r0, r0, lsr #2\n\t"
    "eor     r0, r0, r0, lsr #1\n\t"
    "and     r0, r0, #0x1\n\t"
#endif // USE_PARITY_TABLE

#endif // USE_PARITYSI2

    // Branch to .ERROR if parity fails
    "cbnz    r0, %l[ERROR]\n\t"

    // Mask out parity bits
    "and     r2, r2, #0x00FFFFFF\n\t"

    // Accumulate dot product into result
    "add     r2, %[result], r2, lsl #3\n\t"
    "add     r1, %[vector], r1, lsl #3\n\t"
    "vldr.64 d5, [r2]\n\t"
    "vldr.64 d6, [%[ep], #8]\n\t"
    "vldr.64 d7, [r1]\n\t"
    "vmla.f64        d5, d6, d7\n\t"
    "vstr.64 d5, [r2]\n\t"

    // Increment data pointer, compare to end and branch to loop start
    "add     %[ep], %[ep], #16\n\t"
    "cmp     %[ep], r5\n\t"
    "bne     .LOOP_BODY\n"
    :
    : [result] "r" (result),
      [ep] "r" (matrix.elements),
      [vector] "r" (vector),
      [nnz] "r" (matrix.nnz)
#if USE_PARITY_TABLE
      , [PARITY_TABLE] "r" (PARITY_TABLE)
#endif
    : "cc", "r0", "r1", "r2", "r4", "r5", "d5", "d6", "d7", "memory"
#if USE_PARITYSI2
      , "lr"
#endif
    : ERROR
    );

  return;

ERROR:
  // TODO: Print actual index
  printf("[ECC] error detected at index (something)\n");
  exit(1);
}

static void spmv_sec7(sparse_matrix matrix, double *vector, double *result)
{
  printf("spmv_sec7 not implemented\n");
  exit(1);
}

static void spmv_sec8(sparse_matrix matrix, double *vector, double *result)
{
  printf("spmv_sec8 not implemented\n");
  exit(1);
}

static void spmv_secded(sparse_matrix matrix, double *vector, double *result)
{
  printf("spmv_secded not implemented\n");
  exit(1);
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
