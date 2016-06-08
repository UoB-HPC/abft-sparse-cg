#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>

typedef struct
{
  double value;
  uint32_t column;
} __attribute__((packed)) csr_colval;

typedef struct
{
  unsigned N;
  unsigned nnz;
  uint32_t *cols;
  uint32_t *rows;
  double   *values;
} sparse_matrix;


// Initialize ECC for a sparse matrix
// Defined in the relevant spmv-*.c file
void init_matrix_ecc(sparse_matrix M);

// Sparse matrix-vector product
// Defined in the relevant spmv-*.c file
void spmv(sparse_matrix matrix, double *vector, double *result);

// Flip a specific bit in a matrix element
void flip_bit(csr_colval* element, uint32_t bit);

// Load a sparse matrix from a matrix-market format file
sparse_matrix load_sparse_matrix(const char *matrix_file, int num_blocks);

#endif // COMMON_H
