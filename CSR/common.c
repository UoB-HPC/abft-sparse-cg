#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "ecc.h"

#include "../mmio.h"

typedef struct
{
  uint32_t col;
  uint32_t row;
  double value;
} coo_element;

int compare_matrix_elements(const void *a, const void *b)
{
  coo_element _a = *(coo_element*)a;
  coo_element _b = *(coo_element*)b;

  if (_a.row < _b.row)
  {
    return -1;
  }
  else if (_a.row > _b.row)
  {
    return 1;
  }
  else
  {
    return _a.col - _b.col;
  }
}

sparse_matrix load_sparse_matrix(const char *matrix_file, int num_blocks,
                                 abft_mode mode)
{
  sparse_matrix M = {0};

  FILE *file = fopen(matrix_file, "r");
  if (file == NULL)
  {
    printf("Failed to open '%s'\n", matrix_file);
    exit(1);
  }

  int width, height, nnz;
  mm_read_mtx_crd_size(file, &width, &height, &nnz);
  if (width != height)
  {
    printf("Matrix is not square\n");
    exit(1);
  }

  M.nnz = 0;
  M.rows   = malloc((num_blocks*height+1)*sizeof(uint32_t));
  M.cols   = malloc(num_blocks*2*nnz*sizeof(uint32_t));
  M.values = malloc(num_blocks*2*nnz*sizeof(double));

  coo_element *elements = malloc(2*nnz*sizeof(coo_element));
  int coo_nnz = 0;
  for (int i = 0; i < nnz; i++)
  {
    coo_element element;

    int col, row;

    if (fscanf(file, "%d %d %lg\n", &col, &row, &element.value) != 3)
    {
      printf("Failed to read matrix data\n");
      exit(1);
    }
    col--; /* adjust from 1-based to 0-based */
    row--;

    element.col = col;
    element.row = row;
    elements[coo_nnz] = element;
    coo_nnz++;

    if (element.col == element.row)
      continue;

    element.row = col;
    element.col = row;
    elements[coo_nnz] = element;
    coo_nnz++;
  }

  qsort(elements, coo_nnz, sizeof(coo_element), compare_matrix_elements);

  uint32_t next_row = 0;
  for (int j = 0; j < num_blocks; j++)
  {
    for (int i = 0; i < coo_nnz; i++)
    {
      uint32_t col = elements[i].col + j*width;;
      uint32_t row = elements[i].row + j*height;
      double   val = elements[i].value;

      while (next_row <= row)
      {
        M.rows[next_row++] = M.nnz;
      }

      M.cols[M.nnz]   = col;
      M.values[M.nnz] = val;
      M.nnz++;
    }
  }
  M.rows[height*num_blocks] = M.nnz;

  free(elements);

  M.N = width*num_blocks;
  M.mode = mode;

  // Initialize ECC bits
  for (int i = 0; i < M.nnz; i++)
  {
    csr_colval colval;
    colval.column = M.cols[i];
    colval.value = M.values[i];

    switch (mode)
    {
    case NONE:
    case CONSTRAINTS:
      break;
    case SED:
      colval.column |= ecc_compute_overall_parity(colval) << 31;
      break;
    case SEC7:
      colval.column |= ecc_compute_col8(colval);
      break;
    case SEC8:
    case SECDED:
      colval.column |= ecc_compute_col8(colval);
      colval.column |= ecc_compute_overall_parity(colval) << 24;
      break;
    }

    M.cols[i] = colval.column;
    M.values[i] = colval.value;
  }

  return M;
}

void flip_bit(csr_colval *element, uint32_t bit)
{
  ((uint32_t*)element)[bit/32] ^= 0x1 << (bit % 32);
}
