#include <stdio.h>
#include <stdlib.h>

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

    // Compute overall parity bit for whole codeword
    colval.column |= ecc_compute_overall_parity(colval) << 31;

    M.cols[i] = colval.column;
  }
}

// Sparse matrix vector product
// Multiplies `matrix` by `vector` and stores answer in `result`
void spmv(sparse_matrix matrix, double *vector, double *result)
{
  double zero = 0.0;
  for (unsigned row = 0; row < matrix.N; row++)
  {
    asm goto(
      "add          r5, %[rowptr], %[row], lsl #2\n\t"
      "ldr          r6, [r5]\n\t"
      "ldr          r5, [r5, #4]\n\t"

      "vmov.f64     d5, %P[zero]\n\t"

      ".LOOP_BEGIN:\n\t"
      "cmp          r5, r6\n\t"
      "beq          .LOOP_END\n\t"

      "add     r2, %[colptr], r6, lsl #2\n\t"
      "ldr     r2, [r2]\n\t"

      "add     r1, %[valptr], r6, lsl #3\n\t"
      "vldr.64 d6, [r1]\n\t"

      "vmov    r0, r1, d6\n\t"
      "eor     r0, r0, r1\n\t"
      "eor     r0, r0, r2\n\t"

      "eor     r0, r0, r0, lsr #16\n\t"
      "eor     r0, r0, r0, lsr #8\n\t"
      "eor     r0, r0, r0, lsr #4\n\t"
      "eor     r0, r0, r0, lsr #2\n\t"
      "eor     r0, r0, r0, lsr #1\n\t"
      "and     r0, r0, #0x1\n\t"

      "cbnz    r0, %l[ERROR]\n\t"

      "and     r2, r2, #0x00FFFFFF\n\t"

      "add     r1, %[vecptr], r2, lsl #3\n\t"
      "vldr.64 d7, [r1]\n\t"
      "vmla.f64        d5, d6, d7\n\t"

      "add     r6, r6, #1\n\t"
      "b       .LOOP_BEGIN\n\t"
      ".LOOP_END:\n\t"

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
      : "cc", "r0", "r1", "r2", "r5", "r6", "d5", "d6", "d7", "memory"
      : ERROR
      );
  }

  return;

ERROR:
  printf("[ECC] error detected at index (something)\n");
  exit(1);
}
