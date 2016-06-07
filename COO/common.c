#include <stdio.h>
#include <stdlib.h>

#include "common.h"

#include "../mmio.h"

int compare_matrix_elements(const void *a, const void *b)
{
  matrix_entry _a = *(matrix_entry*)a;
  matrix_entry _b = *(matrix_entry*)b;

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

sparse_matrix load_sparse_matrix(const char *matrix_file, int num_blocks)
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
  M.elements = malloc(num_blocks*2*nnz*sizeof(matrix_entry));
  for (int i = 0; i < nnz; i++)
  {
    matrix_entry element;

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
    M.elements[M.nnz] = element;
    M.nnz++;

    if (element.col == element.row)
      continue;

    element.row = col;
    element.col = row;
    M.elements[M.nnz] = element;
    M.nnz++;
  }

  qsort(M.elements, M.nnz, sizeof(matrix_entry), compare_matrix_elements);

  nnz = M.nnz;
  for (int j = 1; j < num_blocks; j++)
  {
    for (int i = 0; i < nnz; i++)
    {
      matrix_entry element = M.elements[i];
      element.col = element.col + j*width;
      element.row = element.row + j*height;
      M.elements[M.nnz] = element;
      M.nnz++;
    }
  }

  M.N = width*num_blocks;

  return M;
}
