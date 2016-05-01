#include <stdint.h>

// 128-bit matrix element
// Bits  0 to  31 are the colum index
// Bits 32 to  63 are the row index
// Bits 64 to 127 are the floating point value
typedef struct
{
  uint32_t col;
  uint32_t row;
  double value;
} matrix_entry;

typedef struct
{
  unsigned nnz;
  matrix_entry *elements;
} sparse_matrix;


// Initialize ECC for a sparse matrix
// Defined in the relevant spmv-*.c file
void init_matrix_ecc(sparse_matrix M);

// Sparse matrix-vector product
// Defined in the relevant spmv-*.c file
void spmv(sparse_matrix matrix, double *vector, double *result, unsigned N);


// This function will generate/check the 7 parity bits for the given matrix
// element, with the parity bits stored in the high order bits of the column
// index.
//
// This will return a 32-bit integer where the high 7 bits are the generated
// parity bits.
//
// To check a matrix element for errors, simply use this function again, and
// the returned value will be the error 'syndrome' which will be non-zero if
// an error occured.
uint32_t ecc_compute_col8(matrix_entry element);

// This function will correct a single bit-flip in the provided matrix element,
// using the error 'syndrome' generated from a 7-bit parity check.
void     ecc_correct_col8(matrix_entry *element, uint32_t syndrome);

// Compute the overall parity of a 128-bit matrix element
uint32_t ecc_compute_overall_parity(matrix_entry element);
