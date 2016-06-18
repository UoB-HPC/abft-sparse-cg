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
#endif // USE_PARITY_TABLE

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
  double zero = 0.0;
  for (unsigned row = 0; row < matrix.N; row++)
  {
    asm goto(
      // Load row start and row end indices
      "add          r5, %[rowptr], %[row], lsl #2\n\t"
      "ldr          r6, [r5]\n\t"
      "ldr          r5, [r5, #4]\n\t"

      // Zero accumulator
      "vmov.f64     d5, %P[zero]\n\t"

      ".LOOP_BEGIN:\n\t"
      "cmp          r5, r6\n\t"
      "beq          .LOOP_END\n\t"

      // Load column index
      "add     r2, %[colptr], r6, lsl #2\n\t"
      "ldr     r2, [r2]\n\t"

      // Load matrix value
      "add     r1, %[valptr], r6, lsl #3\n\t"
      "vldr.64 d6, [r1]\n\t"

      // *** Parity check starts ***
      // Reduce data to 32-bits in r0
      "ldrd    r0, r1, [r1]\n\t"
      "eor     r0, r0, r1\n\t"
      "eor     r0, r0, r2\n\t"

#if USE_PARITYSI2
      // Use builtin to compute final parity
      "bl      __paritysi2\n\t"
#else
      // Compute final parity manually
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
      "add     r1, %[vecptr], r2, lsl #3\n\t"
      "vldr.64 d7, [r1]\n\t"
      "vmla.f64        d5, d6, d7\n\t"

      // Increment data pointer and branch to loop start
      "add     r6, r6, #1\n\t"
      "b       .LOOP_BEGIN\n\t"
      ".LOOP_END:\n\t"

      // Store accumulator to result vector
      "add     r2, %[resptr], %[row], lsl #3\n\t"
      "vstr.64 d5, [r2]\n\t"

      :
      : [row] "r" (row),
        [rowptr] "r" (matrix.rows),
        [colptr] "r" (matrix.cols),
        [valptr] "r" (matrix.values),
        [vecptr] "r" (vector),
        [resptr] "r" (result),
        [zero] "w" (zero)
#if USE_PARITY_TABLE
        , [PARITY_TABLE] "r" (PARITY_TABLE)
#endif
      : "cc", "r0", "r1", "r2", "r5", "r6", "d5", "d6", "d7", "memory"
#if USE_PARTIYSI2
        , "lr"
#endif
      : ERROR
      );
  }

  return;

ERROR:
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