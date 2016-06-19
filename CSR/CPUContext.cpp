#include "CGContext.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ecc.h"

struct cg_vector
{
  int N;
  double *data;
};

struct cg_matrix
{
  unsigned N;
  unsigned nnz;
  uint32_t *cols;
  uint32_t *rows;
  double   *values;
};

class CPUContext : public CGContext
{
  virtual void generate_ecc_bits(csr_element& element)
  {
  }

  cg_matrix* create_matrix(const uint32_t *columns,
                           const uint32_t *rows,
                           const double *values,
                           int N, int nnz)
  {
    cg_matrix *M = new cg_matrix;

    M->N      = N;
    M->nnz    = nnz;
    M->cols   = new uint32_t[nnz];
    M->rows   = new uint32_t[N+1];
    M->values = new double[nnz];

    uint32_t next_row = 0;
    for (int i = 0; i < nnz; i++)
    {
      csr_element element;
      element.column = columns[i];
      element.value  = values[i];

      generate_ecc_bits(element);

      M->cols[i]   = element.column;
      M->values[i] = element.value;

      while (next_row <= rows[i])
      {
        M->rows[next_row++] = i;
      }
    }
    M->rows[N] = nnz;

    return M;
  }

  void destroy_matrix(cg_matrix *mat)
  {
    delete[] mat->cols;
    delete[] mat->rows;
    delete[] mat->values;
    delete mat;
  }

  cg_vector* create_vector(int N)
  {
    cg_vector *result = new cg_vector;
    result->N    = N;
    result->data = new double[N];
    return result;
  }

  void destroy_vector(cg_vector *vec)
  {
    delete[] vec->data;
    delete vec;
  }

  double* map_vector(cg_vector *v)
  {
    return v->data;
  }

  void unmap_vector(cg_vector *v, double *h)
  {
  }

  void copy_vector(cg_vector *dst, cg_vector *src)
  {
    memcpy(dst->data, src->data, dst->N*sizeof(double));
  }

  double dot(cg_vector *a, cg_vector *b)
  {
    double ret = 0.0;
    for (int i = 0; i < a->N; i++)
    {
      ret += a->data[i] * b->data[i];
    }
    return ret;
  }

  double calc_xr(cg_vector *x, cg_vector *r,
                 cg_vector *p, cg_vector *w,
                 double alpha)
  {
    double ret = 0.0;
    for (int i = 0; i < x->N; i++)
    {
      x->data[i] += alpha * p->data[i];
      r->data[i] -= alpha * w->data[i];

      ret += r->data[i] * r->data[i];
    }
    return ret;
  }

  void calc_p(cg_vector *p, cg_vector *r, double beta)
  {
    for (int i = 0; i < p->N; i++)
    {
      p->data[i] = r->data[i] + beta*p->data[i];
    }
  }

  void spmv(cg_matrix *mat, cg_vector *vec, cg_vector *result)
  {
    for (int row = 0; row < mat->N; row++)
    {
      double tmp = 0.0;

      uint32_t start = mat->rows[row];
      uint32_t end   = mat->rows[row+1];
      for (uint32_t i = start; i < end; i++)
      {
        uint32_t col = mat->cols[i];
        tmp += mat->values[i] * vec->data[col];
      }

      result->data[row] = tmp;
    }
  }

  void inject_bitflip(cg_matrix *mat, BitFlipKind kind, int num_flips)
  {
    int index = rand() % mat->nnz;

    int start = 0;
    int end   = 96;
    if (kind == VALUE)
      end = 64;
    else if (kind == INDEX)
      start = 64;

    for (int i = 0; i < num_flips; i++)
    {
      int bit = (rand() % (end-start)) + start;
      if (bit < 64)
      {
        printf("*** flipping bit %d of value at index %d ***\n", bit, index);
        *((uint64_t*)mat->values+index) ^= 0x1 << (bit % 32);
      }
      else
      {
        bit = bit - 64;
        printf("*** flipping bit %d of column at index %d ***\n", bit, index);
        mat->cols[index] ^= 0x1 << (bit % 32);
      }
    }
  }
};

class CPUContext_Constraints : public CPUContext
{
  void spmv(cg_matrix *mat, cg_vector *vec, cg_vector *result)
  {
    for (int row = 0; row < mat->N; row++)
    {
      double tmp = 0.0;

      uint32_t start = mat->rows[row];
      uint32_t end   = mat->rows[row+1];

      if (end > mat->nnz)
      {
        printf("row size constraint violated for row %d\n", row);
        exit(1);
      }
      if (end < start)
      {
        printf("row order constraint violated for row%d\n", row);
        exit(1);
      }

      for (uint32_t i = start; i < end; i++)
      {
        uint32_t col = mat->cols[i];

        if (col >= mat->N)
        {
          printf("column size constraint violated at index %d\n", i);
          exit(1);
        }
        if (i < end-1)
        {
          if (mat->cols[i+1] <= col)
          {
            printf("column order constraint violated at index %d\n", i);
            exit(1);
          }
        }

        tmp += mat->values[i] * vec->data[col];
      }

      result->data[row] = tmp;
    }
  }
};

class CPUContext_SED : public CPUContext
{
  virtual void generate_ecc_bits(csr_element& element)
  {
    element.column |= ecc_compute_overall_parity(element) << 31;
  }

  void spmv(cg_matrix *mat, cg_vector *vec, cg_vector *result)
  {
    for (int row = 0; row < mat->N; row++)
    {
      double tmp = 0.0;

      uint32_t start = mat->rows[row];
      uint32_t end   = mat->rows[row+1];
      for (uint32_t i = start; i < end; i++)
      {
        csr_element element;
        element.value  = mat->values[i];
        element.column = mat->cols[i];

        // Check overall parity bit
        if (ecc_compute_overall_parity(element))
        {
          printf("[ECC] error detected at index %d\n", i);
          exit(1);
        }

        // Mask out ECC from high order column bits
        element.column &= 0x00FFFFFF;

        tmp += mat->values[i] * vec->data[element.column];
      }

      result->data[row] = tmp;
    }
  }
};

class CPUContext_SEC7 : public CPUContext
{
  virtual void generate_ecc_bits(csr_element& element)
  {
    element.column |= ecc_compute_col8(element);
  }

  void spmv(cg_matrix *mat, cg_vector *vec, cg_vector *result)
  {
    for (int row = 0; row < mat->N; row++)
    {
      double tmp = 0.0;

      uint32_t start = mat->rows[row];
      uint32_t end   = mat->rows[row+1];
      for (uint32_t i = start; i < end; i++)
      {
        csr_element element;
        element.value  = mat->values[i];
        element.column = mat->cols[i];

        // Check ECC
        uint32_t syndrome = ecc_compute_col8(element);
        if (syndrome)
        {
          // Unflip bit
          uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
          ((uint32_t*)(&element))[bit/32] ^= 0x1 << (bit % 32);
          mat->cols[i] = element.column;
          mat->values[i] = element.value;

          printf("[ECC] corrected bit %u at index %d\n", bit, i);
        }

        // Mask out ECC from high order column bits
        element.column &= 0x00FFFFFF;

        tmp += mat->values[i] * vec->data[element.column];
      }

      result->data[row] = tmp;
    }
  }
};

class CPUContext_SEC8 : public CPUContext
{
  virtual void generate_ecc_bits(csr_element& element)
  {
    element.column |= ecc_compute_col8(element);
    element.column |= ecc_compute_overall_parity(element) << 24;
  }

  void spmv(cg_matrix *mat, cg_vector *vec, cg_vector *result)
  {
    for (int row = 0; row < mat->N; row++)
    {
      double tmp = 0.0;

      uint32_t start = mat->rows[row];
      uint32_t end   = mat->rows[row+1];
      for (uint32_t i = start; i < end; i++)
      {
        csr_element element;
        element.value  = mat->values[i];
        element.column = mat->cols[i];

        // Check overall parity bit
        if (ecc_compute_overall_parity(element))
        {
          // Compute error syndrome from hamming bits
          uint32_t syndrome = ecc_compute_col8(element);
          if (syndrome)
          {
            // Unflip bit
            uint32_t bit = ecc_get_flipped_bit_col8(syndrome);
            ((uint32_t*)(&element))[bit/32] ^= 0x1 << (bit % 32);

            printf("[ECC] corrected bit %u at index %d\n", bit, i);
          }
          else
          {
            // Correct overall parity bit
            element.column ^= 0x1 << 24;

            printf("[ECC] corrected overall parity bit at index %d\n", i);
          }
          mat->cols[i] = element.column;
          mat->values[i] = element.value;
        }

        // Mask out ECC from high order column bits
        element.column &= 0x00FFFFFF;

        tmp += mat->values[i] * vec->data[element.column];
      }

      result->data[row] = tmp;
    }
  }
};

namespace
{
  static CGContext::Register<CPUContext> A("cpu", "none");
  static CGContext::Register<CPUContext_Constraints> B("cpu", "constraints");
  static CGContext::Register<CPUContext_SED> C("cpu", "sed");
  static CGContext::Register<CPUContext_SEC7> D("cpu", "sec7");
  static CGContext::Register<CPUContext_SEC8> E("cpu", "sec8");
}
