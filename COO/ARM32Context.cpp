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
    // Initialize result vector to zero
    for (unsigned i = 0; i < mat->N; i++)
      result->data[i] = 0.0;

    coo_element *elements = mat->elements;
    asm volatile(
      // Compute pointer to end of data
      "add     r4, %[elements], %[nnz], lsl #4\n"

      ".LOOP_BODY:\n\t"

      // Load column to r2 and row to r1
      "ldr     r2, [%[elements]]\n\t"
      "ldr     r1, [%[elements], #4]\n\t"

      // *** Parity check starts ***
      // Reduce data to 32-bits in r0
      "eor     r0, r2, r1\n\t"
      "ldr     r3, [%[elements], #8]\n\t"
      "eor     r0, r0, r3\n\t"
      "ldr     r3, [%[elements], #12]\n\t"
      "eor     r0, r0, r3\n\t"

      "eor     r0, r0, r0, lsr #16\n\t"
      "eor     r0, r0, r0, lsr #8\n\t"

      // Lookup final parity from table
      "and     r0, r0, #0xFF\n\t"
      "ldrB    r0, [%[PARITY_TABLE], r0]\n\t"

      // Exit loop if parity fails
      "cbnz    r0, .LOOP_END\n\t"
      // *** Parity check ends ***

      // Mask out parity bits
      "and      r2, r2, #0x00FFFFFF\n\t"

      // Accumulate dot product into result
      "add      r2, %[result], r2, lsl #3\n\t"
      "add      r1, %[vector], r1, lsl #3\n\t"
      "vldr.64  d5, [r2]\n\t"
      "vldr.64  d6, [%[elements], #8]\n\t"
      "vldr.64  d7, [r1]\n\t"
      "vmla.f64 d5, d6, d7\n\t"
      "vstr.64  d5, [r2]\n\t"

      // Increment data pointer, compare to end and branch to loop start
      "add     %[elements], %[elements], #16\n\t"
      "cmp     %[elements], r4\n\t"
      "bne     .LOOP_BODY\n"
      ".LOOP_END:\n"
      : [elements] "+r" (elements)
      : [result] "r" (result->data),
        [vector] "r" (vec->data),
        [nnz] "r" (mat->nnz),
        [PARITY_TABLE] "r" (PARITY_TABLE)
      : "cc", "r0", "r1", "r2", "r3", "r4", "d5", "d6", "d7", "memory"
      );

    if (elements < mat->elements+mat->nnz)
    {
      printf("[ECC] error detected at index %d\n", (elements - mat->elements));
      exit(1);
    }
  }
};

namespace
{
  static CGContext::Register<ARM32Context_SED> A("arm32", "sed");
}
