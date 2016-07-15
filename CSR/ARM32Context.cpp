#include "CPUContext.h"

#include <cstdlib>

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

class ARM32Context_SED : public CPUContext_SED
{
  virtual void spmv(const cg_matrix *mat, const cg_vector *vec, cg_vector *result)
  {
#pragma omp parallel for
    for (unsigned row = 0; row < mat->N; row++)
    {
      double tmp = 0.0;

      int32_t err_index = -1;
      asm (
        // Compute pointers to start of column/value data for this row
        "ldr     r4, [%[rowptr]]\n\t"
        "add     r1, %[cols], r4, lsl #2\n\t"
        "add     r4, %[values], r4, lsl #3\n\t"

        // Compute pointer to end of column data for this row
        "ldr     r5, [%[rowptr], #4]\n\t"
        "add     r5, %[cols], r5, lsl #2\n\t"

        ".LOOP_BODY:\n\t"
        // Check if we've reached the end of this row
        "cmp     r1, r5\n\t"
        "beq     .LOOP_END\n"

        // *** Parity check starts ***
        // Reduce data to 32-bits in r0
        "ldr     r0, [r4]\n\t"
        "ldr     r2, [r4, #4]\n\t"
        "eor     r0, r0, r2\n\t"

        // Load column into r2
        "ldr     r2, [r1]\n\t"

        // Continue reducing data to 32-bits
        "eor     r0, r0, r2\n\t"
        "eor     r0, r0, r0, lsr #16\n\t"
        "eor     r0, r0, r0, lsr #8\n\t"

        // Lookup final parity from table
        "and     r0, r0, #0xFF\n\t"
        "ldrB    r0, [%[PARITY_TABLE], r0]\n\t"

        // Exit loop if parity fails
        "cbz    r0, .NO_ERROR\n\t"
        "sub    %[err_index], r1, %[cols]\n\t"
        "b      .LOOP_END\n"
        ".NO_ERROR:\n\t"
        // *** Parity check ends ***

        // Mask out parity bits
        "and      r2, r2, #0x00FFFFFF\n\t"

        // Accumulate dot product into result
        "add      r2, %[vector], r2, lsl #3\n\t"
        "vldr.64  d6, [r4]\n\t"
        "vldr.64  d7, [r2]\n\t"
        "vmla.f64 %P[tmp], d6, d7\n\t"

        // Increment data pointer, compare to end and branch to loop start
        "add     r4, #8\n\t"
        "add     r1, #4\n\t"
        "b       .LOOP_BODY\n"
        ".LOOP_END:\n"
        : [tmp] "+w" (tmp), [err_index] "+r" (err_index)
        : [rowptr] "r" (mat->rows+row),
          [cols] "r" (mat->cols),
          [values] "r" (mat->values),
          [vector] "r" (vec->data),
          [PARITY_TABLE] "r" (PARITY_TABLE)
        : "cc", "r0", "r1", "r2", "r4", "r5", "d6", "d7"
        );

      result->data[row] = tmp;

      if (err_index >= 0)
      {
        printf("[ECC] error detected at index %d\n", err_index/4);
        exit(1);
      }
    }
  }
};

namespace
{
  static CGContext::Register<ARM32Context_SED> A("arm32", "sed");
}
